"""
Protein–Ligand Interactions (PLI): receptor/ligand loading, 3D contact detection,
PyMOL .pml generation, structure-based pharmacophores, and PoseView-style 2D diagrams.
"""

import os
from typing import Optional, Union
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from pharmacophore.constants import (
    AMINO_ACID_ALIASES,
    DEFAULT_POCKET_CUTOFF,
    METAL_RESIDUES,
    STANDARD_AMINO_ACIDS,
    WATER_RESIDUES,
)


class PLI:
    def __init__(
        self,
        protein_pdb: str,
        ligand_sdf: Union[str, os.PathLike, Chem.Mol],
        pose_idx: int = 0,
        template_smi: Optional[str] = None,
        pocket_cutoff: float = DEFAULT_POCKET_CUTOFF,
        keep_waters: bool = False,
        keep_metals: bool = True,
        add_hs: bool = False,
    ):
        """
        Initialize a Protein–Ligand Interaction (PLI) analyzer from a receptor .pdb and ligand .sdf (or Chem.Mol).
        :param protein_pdb: str
            File path to the receptor structure in .pdb format.
        :param ligand_sdf: Union[str, os.PathLike, Chem.Mol]
            File path to the 3D ligand .sdf file, or an in-memory RDKit Chem.Mol with 3D coordinates.
        :param pose_idx: int
            Zero-based pose index to load when ligand_sdf is a multi-molecule .sdf file. Defaults to 0.
        :param template_smi: Optional[str]
            Optional SMILES string used to restore bond orders and aromaticity on docked ligand poses while
            preserving 3D coordinates.
        :param pocket_cutoff: float
            Heavy-atom distance cutoff in Ångströms used to extract binding-pocket residues and validate
            coordinate-frame alignment. Defaults to 6.0 Å.
        :param keep_waters: bool
            Whether to retain crystallographic water molecules (HOH/WAT) inside the binding pocket. Defaults to False.
        :param keep_metals: bool
            Whether to retain metal ions (e.g., MG, ZN, FE) inside the binding pocket. Defaults to True.
        :param add_hs: bool
            Optional convenience flag (ProLIF-style). If True, calls RDKit Chem.AddHs(..., addCoords=True) on the
            receptor and ligand after loading. Defaults to False (preserving existing hydrogens as-is).
        """
        if pocket_cutoff <= 0:
            raise ValueError(f"pocket_cutoff must be positive, got {pocket_cutoff}.")

        self.protein_pdb = str(protein_pdb)
        self.ligand_source = ligand_sdf
        self.pose_idx = pose_idx
        self.template_smi = template_smi
        self.pocket_cutoff = float(pocket_cutoff)
        self.keep_waters = bool(keep_waters)
        self.keep_metals = bool(keep_metals)
        self.add_hs = bool(add_hs)

        self.protein = _load_protein(self.protein_pdb, add_hs=self.add_hs)
        self.ligand = _load_ligand(
            ligand_sdf,
            pose_idx=self.pose_idx,
            template_smi=self.template_smi,
            add_hs=self.add_hs,
        )
        (
            self.ligand_no_h,
            self.full_to_heavy_idx,
            self.heavy_to_full_idx,
        ) = _build_heavy_atom_map(self.ligand)

        self.pocket_residues = self._extract_pocket_residues()
        self.interactions_df = None

    def _extract_pocket_residues(self) -> list[dict]:
        """
        Filter receptor atoms (removing waters/unrelated HETATMs according to settings), validate 3D coordinate
        alignment against the ligand, and collect all residues within self.pocket_cutoff Å.
        :return: list[dict]
            Ordered list of pocket residue metadata dictionaries.
        """
        lig_heavy_coords = self.ligand_no_h.GetConformer().GetPositions()
        if len(lig_heavy_coords) == 0:
            raise ValueError("Ligand contains no heavy atoms.")

        prot_conf = self.protein.GetConformer()
        prot_positions = prot_conf.GetPositions()

        residue_atoms: dict[tuple[str, int, str, str], dict] = {}
        kept_heavy_atom_indices: list[int] = []
        kept_heavy_res_keys: list[tuple[str, int, str, str]] = []

        for atom in self.protein.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is None:
                continue

            raw_res_name = info.GetResidueName().strip().upper()
            res_num = int(info.GetResidueNumber())
            ins_code = info.GetInsertionCode().strip()
            chain_id = info.GetChainId().strip()
            is_het = bool(info.GetIsHeteroAtom())

            is_water = raw_res_name in WATER_RESIDUES
            is_metal = raw_res_name in METAL_RESIDUES
            is_amino_acid = (raw_res_name in STANDARD_AMINO_ACIDS) or (not is_het and not is_water and not is_metal)

            if is_water:
                if not self.keep_waters:
                    continue
            elif is_metal:
                if not self.keep_metals:
                    continue
            elif not is_amino_acid:
                # Exclude unrelated HETATM co-solvents, buffers, or co-crystallized ligands
                continue

            res_key = (chain_id, res_num, ins_code, raw_res_name)
            if res_key not in residue_atoms:
                canonical_name = AMINO_ACID_ALIASES.get(raw_res_name, raw_res_name)
                residue_atoms[res_key] = {
                    "chain": chain_id,
                    "residue_number": res_num,
                    "insertion_code": ins_code,
                    "residue_name": canonical_name,
                    "pdb_residue_name": raw_res_name,
                    "is_water": is_water,
                    "is_metal": is_metal,
                    "atom_indices": [],
                    "heavy_atom_indices": [],
                }

            idx = atom.GetIdx()
            residue_atoms[res_key]["atom_indices"].append(idx)
            if atom.GetAtomicNum() > 1:
                residue_atoms[res_key]["heavy_atom_indices"].append(idx)
                kept_heavy_atom_indices.append(idx)
                kept_heavy_res_keys.append(res_key)

        if not kept_heavy_atom_indices:
            raise ValueError(
                f"No valid receptor heavy atoms found in {self.protein_pdb} after filtering."
            )

        prot_heavy_coords = prot_positions[kept_heavy_atom_indices]
        # Pairwise distance matrix: shape (N_prot_heavy, N_lig_heavy)
        dists = np.linalg.norm(
            prot_heavy_coords[:, None, :] - lig_heavy_coords[None, :, :],
            axis=-1,
        )
        min_per_prot_atom = dists.min(axis=1)
        overall_min_dist = float(min_per_prot_atom.min())

        if overall_min_dist >= self.pocket_cutoff:
            raise ValueError(
                f"No receptor atoms found within {self.pocket_cutoff:.1f} Å of the ligand "
                f"(closest heavy-atom distance: {overall_min_dist:.2f} Å). Ensure the protein .pdb "
                f"and ligand .sdf share the same 3D binding-site coordinate frame."
            )

        res_min_dist: dict[tuple[str, int, str, str], float] = {}
        for res_key, d in zip(kept_heavy_res_keys, min_per_prot_atom):
            d_val = float(d)
            if d_val <= self.pocket_cutoff:
                if res_key not in res_min_dist or d_val < res_min_dist[res_key]:
                    res_min_dist[res_key] = d_val

        pocket_list = []
        for res_key in sorted(res_min_dist.keys(), key=lambda k: (k[0], k[1], k[2], k[3])):
            entry = dict(residue_atoms[res_key])
            entry["min_distance"] = round(res_min_dist[res_key], 3)
            pocket_list.append(entry)

        return pocket_list

    def pocket_summary(self) -> pd.DataFrame:
        """
        Return a summary DataFrame of all receptor residues within pocket_cutoff of the ligand.
        :return: pd.DataFrame
            Columns: residue_name, pdb_residue_name, residue_number, chain, min_distance, num_heavy_atoms,
            is_water, is_metal.
        """
        rows = [
            {
                "residue_name": r["residue_name"],
                "pdb_residue_name": r["pdb_residue_name"],
                "residue_number": r["residue_number"],
                "chain": r["chain"],
                "min_distance": r["min_distance"],
                "num_heavy_atoms": len(r["heavy_atom_indices"]),
                "is_water": r["is_water"],
                "is_metal": r["is_metal"],
            }
            for r in self.pocket_residues
        ]
        return pd.DataFrame(
            rows,
            columns=[
                "residue_name",
                "pdb_residue_name",
                "residue_number",
                "chain",
                "min_distance",
                "num_heavy_atoms",
                "is_water",
                "is_metal",
            ],
        )


def _load_protein(protein_pdb: str, add_hs: bool = False) -> Chem.Mol:
    """
    Load a receptor .pdb structure with RDKit while preserving existing hydrogens and handling minor valence quirks.
    """
    if not os.path.isfile(protein_pdb):
        raise FileNotFoundError(f"Protein PDB file not found: {protein_pdb}")
    if os.path.getsize(protein_pdb) == 0:
        raise ValueError(f"Protein PDB file is empty: {protein_pdb}")

    mol = Chem.MolFromPDBFile(protein_pdb, removeHs=False, sanitize=True)
    if mol is None:
        mol = Chem.MolFromPDBFile(protein_pdb, removeHs=False, sanitize=False)
        if mol is None:
            raise ValueError(f"RDKit failed to parse protein PDB file: {protein_pdb}")
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
            catchErrors=True,
        )
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)

    if mol.GetNumConformers() == 0:
        raise ValueError(f"Protein PDB file contains no 3D coordinates: {protein_pdb}")

    if add_hs:
        mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)

    return mol


def _load_ligand(
    ligand_sdf: Union[str, os.PathLike, Chem.Mol],
    pose_idx: int = 0,
    template_smi: Optional[str] = None,
    add_hs: bool = False,
) -> Chem.Mol:
    """
    Load a 3D ligand from an .sdf file (at pose_idx) or copy an in-memory RDKit Chem.Mol.
    """
    if isinstance(ligand_sdf, Chem.Mol):
        mol = Chem.Mol(ligand_sdf)
    elif isinstance(ligand_sdf, (str, os.PathLike)):
        sdf_path = str(ligand_sdf)
        if not os.path.isfile(sdf_path):
            raise FileNotFoundError(f"Ligand SDF file not found: {sdf_path}")
        supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
        num_records = len(supplier)
        if num_records == 0:
            raise ValueError(f"No molecule records found in SDF file: {sdf_path}")
        if pose_idx < 0 or pose_idx >= num_records:
            raise IndexError(
                f"pose_idx={pose_idx} out of range for {sdf_path} ({num_records} record(s))."
            )
        mol = supplier[pose_idx]
        if mol is None:
            raise ValueError(f"RDKit failed to parse ligand at pose_idx={pose_idx} in {sdf_path}.")
    else:
        raise TypeError("ligand_sdf must be a path to an .sdf file or an RDKit Chem.Mol.")

    if mol.GetNumConformers() == 0:
        raise ValueError(
            "Ligand has no 3D conformer. Supply a 3D bound/docked SDF or embed 3D coordinates first."
        )

    if template_smi is not None:
        template = Chem.MolFromSmiles(template_smi)
        if template is None:
            raise ValueError(f"Invalid template_smi SMILES string: {template_smi}")
        had_explicit_hs = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
        mol_no_h = Chem.RemoveHs(mol)
        mol = AllChem.AssignBondOrdersFromTemplate(template, mol_no_h)
        if had_explicit_hs or add_hs:
            mol = Chem.AddHs(mol, addCoords=True)
    elif add_hs:
        mol = Chem.AddHs(mol, addCoords=True)

    return mol


def _build_heavy_atom_map(mol: Chem.Mol) -> tuple[Chem.Mol, dict[int, int], dict[int, int]]:
    """
    Create a heavy-atom-only copy of mol (removing explicit hydrogens) and return bidirectional atom index maps.
    :return: (mol_no_h, full_to_heavy_idx, heavy_to_full_idx)
    """
    tagged = Chem.Mol(mol)
    for atom in tagged.GetAtoms():
        atom.SetIntProp("_pli_full_idx", atom.GetIdx())

    mol_no_h = Chem.RemoveHs(tagged)
    full_to_heavy: dict[int, int] = {}
    heavy_to_full: dict[int, int] = {}

    for atom in mol_no_h.GetAtoms():
        full_idx = atom.GetIntProp("_pli_full_idx")
        heavy_idx = atom.GetIdx()
        full_to_heavy[full_idx] = heavy_idx
        heavy_to_full[heavy_idx] = full_idx
        atom.ClearProp("_pli_full_idx")

    return mol_no_h, full_to_heavy, heavy_to_full
