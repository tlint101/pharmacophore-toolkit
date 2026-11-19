"""
Protein–Ligand Interactions (PLI): receptor/ligand loading, 3D contact detection,
PyMOL .pml generation, structure-based pharmacophores, and protein-ligand interaction 2D diagrams.
"""

import io
import json
import math
import os
import re
from typing import Optional, Union
import contourpy
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from PIL import Image
import py3Dmol
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Geometry import Point3D

from pharmacophore.constants import (
    AMINO_ACID_ALIASES,
    BACKBONE_ATOMS,
    DEFAULT_CATION_PI_ANGLE,
    DEFAULT_CATION_PI_DIST,
    DEFAULT_CLASH_DIST,
    DEFAULT_HBOND_ACCEPTOR_ANGLE,
    DEFAULT_HBOND_ANGLE,
    DEFAULT_HBOND_DIST,
    DEFAULT_HBOND_H_DIST,
    DEFAULT_HBOND_HEAVY_ANGLE,
    DEFAULT_HYDROPHOBIC_DIST,
    DEFAULT_IONIC_DIST,
    DEFAULT_METAL_DIST,
    DEFAULT_PI_EDGE_ANGLE,
    DEFAULT_PI_FACE_ANGLE,
    DEFAULT_PI_FACE_DIST,
    DEFAULT_PI_STACKING_DIST,
    DEFAULT_POCKET_CUTOFF,
    DEFAULT_VDW_TOLERANCE,
    DEFAULT_XBOND_ACCEPTOR_ANGLE,
    DEFAULT_XBOND_DIST,
    DEFAULT_XBOND_DONOR_ANGLE,
    FEATURE_COLORS,
    HYDROPHOBIC_RESIDUES,
    INTERACTIVE_COLORS,
    METAL_ATOMIC_NUMS,
    METAL_RESIDUES,
    PLI_COLORS,
    PML_COLORS,
    PROLIG_COLORS,
    PROLIG_DEFAULT_INTERACTIONS,
    PROLIG_INTERACTION_ORDER,
    PROTEIN_AROMATIC_RESIDUES,
    PROTEIN_AROMATIC_RING_ATOMS,
    PROTEIN_NEG_ATOMS,
    PROTEIN_POS_ATOMS,
    STANDARD_AMINO_ACIDS,
    WATER_RESIDUES,
    color_convert,
    format_pymol_split_spheres_block,
    is_hbond_acceptor_nitrogen,
)

_POS_CHARGED_N_QUERY = Chem.MolFromSmarts("[#7;+;!$([N+]([O-]))]")
_POS_AMINE_QUERY = Chem.MolFromSmarts("[$([NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(Na)])]")
_NEG_ACID_O_QUERY = Chem.MolFromSmarts("[$([OX1-,OX2H1][CX3,SX4,PX4]=O)]")

_SUMMARY_COLUMNS = [
    "residue_name",
    "residue_number",
    "chain",
    "is_backbone",
    "interaction_type",
    "ligand_atom_indices",
    "protein_atom_name",
    "distance",
]

_DETAILED_COLUMNS = _SUMMARY_COLUMNS + [
    "angle",
    "ligand_coords",
    "protein_coords",
]


class PLI:
    def __init__(
        self,
        protein: str,
        ligand: Union[str, os.PathLike, Chem.Mol, list[Chem.Mol], tuple[Chem.Mol, ...]],
        pocket_cutoff: float = DEFAULT_POCKET_CUTOFF,
        keep_waters: bool = False,
        keep_metals: bool = True,
        add_hs: bool = False,
    ):
        """
        Initialize a Protein–Ligand Interaction (PLI) analyzer from a receptor. A protein structure (.pdb) and the ligand
        file (.sdf) should be given. The protein structure is not modified. If waters are included, they are removed.
        Metals are kept by default and non-protein molecules (solvents or additonal co-crystsallized ligands) are
        ignored. For ligands, any existing hdyrogens, charges, and 3D coordinates are kept. No hydrogens are added unless
        specified using "add_hs".
        :param protein: str
            File path to the receptor structure in ``.pdb`` format.
        :param ligand: Union[str, os.PathLike, Chem.Mol, list[Chem.Mol], tuple[Chem.Mol, ...]]
            File path to a single- or multi-molecule 3D ligand ``.sdf`` file, an in-memory RDKit ``Chem.Mol``, or
            a non-empty sequence of 3D ``Chem.Mol`` objects. All loaded ligands are stored and scored by default in
            ``to_df()`` and ``plot_barcode()``.
        :param pocket_cutoff: float
            Heavy-atom distance cutoff in Ångströms used to define the binding pocket. Defaults to ``6.0``.
        :param keep_waters: bool
            Whether to keep crystallographic waters in the binding pocket. Defaults to ``False``.
        :param keep_metals: bool
            Whether to keep metal ions in the binding pocket. Defaults to ``True``.
        :param add_hs: bool
            If ``True``, adds 3D hydrogens to both protein and ligand with RDKit. Defaults to ``False``.
        """
        if pocket_cutoff <= 0:
            raise ValueError(f"pocket_cutoff must be positive, got {pocket_cutoff}.")

        self.protein_pdb = str(protein)
        self.ligand_source = ligand
        self.pocket_cutoff = float(pocket_cutoff)
        self.keep_waters = bool(keep_waters)
        self.keep_metals = bool(keep_metals)
        self.add_hs = bool(add_hs)

        self.protein = _load_protein(self.protein_pdb, add_hs=self.add_hs)
        self._receptor_index = self._index_receptor_residues()

        loaded_mols = _load_ligands(ligand, add_hs=self.add_hs)
        single_ligand = len(loaded_mols) == 1
        entries: list[dict] = []
        for idx, mol_full in enumerate(loaded_mols):
            mol_no_h, f2h, h2f = _build_heavy_atom_map(mol_full)
            pocket_res = self._extract_pocket_residues(
                mol_no_h, validate_alignment=single_ligand
            )
            prop_name = mol_full.GetProp("_Name").strip() if mol_full.HasProp("_Name") else ""
            entries.append(
                {
                    "name": prop_name or f"Mol_{idx}",
                    "ligand": mol_full,
                    "ligand_no_h": mol_no_h,
                    "full_to_heavy_idx": f2h,
                    "heavy_to_full_idx": h2f,
                    "pocket_residues": pocket_res,
                }
            )

        if not single_ligand and not any(e["pocket_residues"] for e in entries):
            self._extract_pocket_residues(entries[0]["ligand_no_h"], validate_alignment=True)

        self._ligand_entries = entries
        self.ligands = [e["ligand"] for e in entries]
        self.mol_names = [e["name"] for e in entries]
        self._active_idx = -1
        self._select_ligand(0)

    def _resolve_ligand_index(self, mol: Union[int, str] = 0) -> int:
        """Resolve a 0-based index or molecule name into a validated ligand entry index."""
        n_ligs = len(self._ligand_entries)
        if isinstance(mol, int) and not isinstance(mol, bool):
            idx = mol if mol >= 0 else n_ligs + mol
            if idx < 0 or idx >= n_ligs:
                raise IndexError(f"mol={mol} out of range for {n_ligs} loaded ligand(s).")
            return idx
        if isinstance(mol, str):
            if mol in self.mol_names:
                return self.mol_names.index(mol)
            lower_map = {name.lower(): i for i, name in enumerate(self.mol_names)}
            if mol.lower() in lower_map:
                return lower_map[mol.lower()]
            raise KeyError(f"Ligand {mol!r} not found. Available ligands: {self.mol_names}")
        raise TypeError("mol must be a 0-based integer index or a molecule name string.")

    def _select_ligand(self, mol: Union[int, str] = 0) -> int:
        """Bind active single-ligand attributes and reset cached contact state when switching ligands."""
        idx = self._resolve_ligand_index(mol)
        if idx != self._active_idx:
            entry = self._ligand_entries[idx]
            self.ligand = entry["ligand"]
            self.ligand_no_h = entry["ligand_no_h"]
            self.full_to_heavy_idx = entry["full_to_heavy_idx"]
            self.heavy_to_full_idx = entry["heavy_to_full_idx"]
            self.pocket_residues = entry["pocket_residues"]
            self._active_idx = idx
            self._interaction_records = []
            self.interactions_df = None
            self.features_list = None
        return idx

    def _index_receptor_residues(self) -> dict:
        """
        Parse and cache receptor residue metadata and heavy-atom 3D coordinates once per PLI instance.
        :return: dict
            Dictionary with 'residue_atoms', 'kept_heavy_atom_indices', 'kept_heavy_res_keys',
            and 'prot_heavy_coords'.
        """
        prot_positions = self.protein.GetConformer().GetPositions()
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

        return {
            "residue_atoms": residue_atoms,
            "kept_heavy_atom_indices": kept_heavy_atom_indices,
            "kept_heavy_res_keys": kept_heavy_res_keys,
            "prot_heavy_coords": prot_positions[kept_heavy_atom_indices],
        }

    def _extract_pocket_residues(
        self,
        lig_no_h: Optional[Chem.Mol] = None,
        validate_alignment: bool = True,
    ) -> list[dict]:
        """
        Validate 3D coordinate alignment against a heavy-atom ligand and collect all cached receptor residues
        within self.pocket_cutoff Å.
        :param lig_no_h: Optional[Chem.Mol]
            Heavy-atom ligand molecule. Defaults to self.ligand_no_h.
        :param validate_alignment: bool
            If True (default), raise ValueError when no receptor atom is within pocket_cutoff Å.
        :return: list[dict]
            Ordered list of pocket residue metadata dictionaries.
        """
        target_lig = self.ligand_no_h if lig_no_h is None else lig_no_h
        lig_heavy_coords = target_lig.GetConformer().GetPositions()
        if len(lig_heavy_coords) == 0:
            raise ValueError("Ligand contains no heavy atoms.")

        idx_cache = self._receptor_index
        prot_heavy_coords = idx_cache["prot_heavy_coords"]
        # Pairwise distance matrix: shape (N_prot_heavy, N_lig_heavy)
        dists = np.linalg.norm(
            prot_heavy_coords[:, None, :] - lig_heavy_coords[None, :, :],
            axis=-1,
        )
        min_per_prot_atom = dists.min(axis=1)
        overall_min_dist = float(min_per_prot_atom.min())

        if overall_min_dist >= self.pocket_cutoff:
            if not validate_alignment:
                return []
            raise ValueError(
                f"No receptor atoms found within {self.pocket_cutoff:.1f} Å of the ligand "
                f"(closest heavy-atom distance: {overall_min_dist:.2f} Å). Ensure the protein .pdb "
                f"and ligand .sdf share the same 3D binding-site coordinate frame."
            )

        res_min_dist: dict[tuple[str, int, str, str], float] = {}
        for res_key, d in zip(idx_cache["kept_heavy_res_keys"], min_per_prot_atom):
            d_val = float(d)
            if d_val <= self.pocket_cutoff:
                if res_key not in res_min_dist or d_val < res_min_dist[res_key]:
                    res_min_dist[res_key] = d_val

        residue_atoms = idx_cache["residue_atoms"]
        pocket_list = []
        for res_key in sorted(res_min_dist.keys(), key=lambda k: (k[0], k[1], k[2], k[3])):
            entry = dict(residue_atoms[res_key])
            entry["min_distance"] = round(res_min_dist[res_key], 3)
            pocket_list.append(entry)

        return pocket_list

    def pocket_summary(self, mol: Union[int, str] = 0) -> pd.DataFrame:
        """
        Return a summary DataFrame of all receptor residues within pocket_cutoff of the selected ligand.
        :param mol: Union[int, str]
            0-based index or molecule name specifying which loaded ligand to summarize. Defaults to ``0``.
        :return: pd.DataFrame
            Columns: residue_name, pdb_residue_name, residue_number, chain, min_distance, num_heavy_atoms,
            is_water, is_metal.
        """
        entry = self._ligand_entries[self._resolve_ligand_index(mol)]
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
            for r in entry["pocket_residues"]
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

    def _detect(
        self,
        include_hbonds: bool = True,
        include_hydrophobic: bool = True,
        include_aromatic: bool = True,
        include_ionic: bool = True,
        include_clashes: bool = False,
        detailed: bool = False,
        hbond_dist: float = DEFAULT_HBOND_DIST,
        hbond_angle: float = DEFAULT_HBOND_ANGLE,
        hydrophobic_dist: float = DEFAULT_HYDROPHOBIC_DIST,
        ionic_dist: float = DEFAULT_IONIC_DIST,
        pi_stacking_dist: float = DEFAULT_PI_STACKING_DIST,
        cation_pi_dist: float = DEFAULT_CATION_PI_DIST,
        clash_dist: float = DEFAULT_CLASH_DIST,
    ) -> pd.DataFrame:
        """
        Detect 3D protein–ligand interactions within the extracted binding pocket.
        :param include_hbonds: bool
            Detect hydrogen bonds (HBond_Donor / HBond_Acceptor from the ligand's perspective).
        :param include_hydrophobic: bool
            Detect hydrophobic contacts grouped per contacting residue (excluding residues with specific
            polar/aromatic contacts).
        :param include_aromatic: bool
            Detect aromatic pi-pi stacking and cation-pi interactions.
        :param include_ionic: bool
            Detect electrostatic salt bridges / ionic interactions.
        :param include_clashes: bool
            Detect non-H-bonding heavy-atom steric clashes below clash_dist Å.
        :param detailed: bool
            If True, also include angle, ligand_coords, and protein_coords.
        :param hbond_dist: float
            Maximum heavy-atom donor–acceptor distance (Å) for hydrogen bonds.
        :param hbond_angle: float
            Minimum D-H...A angle (degrees) when explicit donor hydrogen are present.
        :param hydrophobic_dist: float
            Maximum heavy-atom distance (Å) for hydrophobic contacts.
        :param ionic_dist: float
            Maximum heavy-atom distance (Å) for salt bridges / ionic interactions.
        :param pi_stacking_dist: float
            Maximum ring-centroid distance (Å) for aromatic pi-pi stacking.
        :param cation_pi_dist: float
            Maximum distance (Å) between a cation and an aromatic ring centroid for cation-pi contacts.
        :param clash_dist: float
            Distance threshold (Å) below which non-H-bonding heavy-atom contacts count as steric clashes.
        :return: pd.DataFrame
            Detected interactions sorted by chain, residue number, and interaction type.
        """
        lig_feats = self._classify_ligand_features()
        rows: list[dict] = []
        hbond_pairs: set[tuple[tuple, int]] = set()
        specific_contact_residues: set[tuple] = set()

        if include_hbonds:
            hb_rows, hbond_pairs = self._detect_hbonds(lig_feats, hbond_dist, hbond_angle)
            rows.extend(hb_rows)
            specific_contact_residues.update(r["_res_key"] for r in hb_rows)

        if include_ionic:
            ionic_rows = self._detect_ionic(lig_feats, ionic_dist, hbond_pairs)
            rows.extend(ionic_rows)
            specific_contact_residues.update(r["_res_key"] for r in ionic_rows)

        if include_aromatic:
            arom_rows = self._detect_aromatic(lig_feats, pi_stacking_dist, cation_pi_dist)
            rows.extend(arom_rows)
            specific_contact_residues.update(r["_res_key"] for r in arom_rows)

        if include_hydrophobic:
            hydro_rows = self._detect_hydrophobic(
                lig_feats, hydrophobic_dist, specific_contact_residues
            )
            rows.extend(hydro_rows)

        if include_clashes:
            clash_rows = self._detect_clashes(clash_dist, hbond_pairs)
            rows.extend(clash_rows)

        for row in rows:
            row.pop("_res_key", None)

        rows.sort(key=lambda r: (r["chain"], r["residue_number"], r["interaction_type"]))
        self._interaction_records = rows

        columns = _DETAILED_COLUMNS if detailed else _SUMMARY_COLUMNS
        df = pd.DataFrame(rows, columns=columns)
        self.interactions_df = df
        return df

    def to_df(
        self,
        mol_names: Optional[list[str]] = None,
        interactions: Union[str, list[str]] = "default",
        count: Union[str, bool] = "ligand",
        show_chain: Union[bool, str] = "auto",
        drop_empty: bool = True,
        **cutoffs,
    ) -> pd.DataFrame:
        """
        Generate a protein–ligand interaction count table across all loaded ligands.
        :param mol_names: Optional[list[str]]
            Optional display names corresponding to each loaded ligand row. When None, uses each molecule's
            ``_Name`` property or falls back to ``'Mol_0'``, ``'Mol_1'``, ...
        :param interactions: Union[str, list[str]]
            Interaction types to compute: 'default' (HBDonor, HBAcceptor, PiStacking, Hydrophobic, VdWContact),
            'all' (all supported classes with PiStacking combined), or an explicit list of interaction
            names from PROLIG_INTERACTION_ORDER (e.g., ['HBDonor', 'HBAcceptor', 'FaceToFace', 'EdgeToFace']).
        :param count: Union[str, bool]
            Contact counting mode: 'ligand' (default, unique contacting ligand atoms or rings per residue),
            'pairs' (ProLIG count=True total atom/ring pairs), or 'binary' (0/1 presence). True/False map to
            'ligand'/'binary'.
        :param show_chain: Union[bool, str]
            Whether to append chain IDs to residue column headers ('auto' appends chains only when multiple
            chains appear in the pocket; True always appends; False omits).
        :param drop_empty: bool
            If True (default), omit interaction columns whose count is 0 across all rows.
        :param cutoffs: dict
            Optional geometric cutoff overrides: hbond_dist (3.5), hbond_angle (120.0), xbond_dist (3.5),
            xbond_donor_angle (120.0), xbond_acceptor_angle (90.0), hydrophobic_dist (4.0), ionic_dist (5.0),
            pi_stacking_dist (5.5), pi_face_dist (5.0), cation_pi_dist (6.0), cation_pi_angle (30.0),
            metal_dist (2.8), vdw_tolerance (0.0).
        :return: pd.DataFrame
            DataFrame with 'mol_name', 'smiles', and grouped <Interaction>_<RES><NUM> integer count columns.
        """
        if count is True:
            count_mode = "ligand"
        elif count is False:
            count_mode = "binary"
        elif isinstance(count, str) and count.lower() in ("ligand", "pairs", "binary"):
            count_mode = count.lower()
        else:
            raise ValueError(
                f"Invalid count={count!r}. Expected 'ligand', 'pairs', 'binary', True, or False."
            )

        if show_chain not in (True, False, "auto"):
            raise ValueError(f"Invalid show_chain={show_chain!r}. Expected True, False, or 'auto'.")

        if isinstance(interactions, str):
            preset = interactions.strip().lower()
            if preset == "default":
                active = list(PROLIG_DEFAULT_INTERACTIONS)
            elif preset == "all":
                active = [
                    name
                    for name in PROLIG_INTERACTION_ORDER
                    if name not in ("FaceToFace", "EdgeToFace")
                ]
            else:
                raise ValueError(
                    f"Unknown interactions preset {interactions!r}. Use 'default', 'all', or a list from "
                    f"{list(PROLIG_INTERACTION_ORDER)}."
                )
        else:
            req = list(interactions)
            unknown = [name for name in req if name not in PROLIG_INTERACTION_ORDER]
            if unknown:
                raise ValueError(
                    f"Unknown interaction type(s): {unknown}. Supported: {list(PROLIG_INTERACTION_ORDER)}."
                )
            req_set = set(req)
            active = [name for name in PROLIG_INTERACTION_ORDER if name in req_set]

        default_cutoffs = {
            "hbond_dist": float(DEFAULT_HBOND_DIST),
            "hbond_angle": float(DEFAULT_HBOND_ANGLE),
            "xbond_dist": float(DEFAULT_XBOND_DIST),
            "xbond_donor_angle": float(DEFAULT_XBOND_DONOR_ANGLE),
            "xbond_acceptor_angle": float(DEFAULT_XBOND_ACCEPTOR_ANGLE),
            "hydrophobic_dist": float(DEFAULT_HYDROPHOBIC_DIST),
            "ionic_dist": float(DEFAULT_IONIC_DIST),
            "pi_stacking_dist": float(DEFAULT_PI_STACKING_DIST),
            "pi_face_dist": float(DEFAULT_PI_FACE_DIST),
            "cation_pi_dist": float(DEFAULT_CATION_PI_DIST),
            "cation_pi_angle": float(DEFAULT_CATION_PI_ANGLE),
            "metal_dist": float(DEFAULT_METAL_DIST),
            "vdw_tolerance": float(DEFAULT_VDW_TOLERANCE),
        }
        unknown_cutoffs = sorted(set(cutoffs.keys()) - set(default_cutoffs.keys()))
        if unknown_cutoffs:
            raise ValueError(
                f"Unknown cutoff argument(s): {unknown_cutoffs}. Supported: {sorted(default_cutoffs.keys())}."
            )
        resolved_cutoffs = dict(default_cutoffs)
        for key, val in cutoffs.items():
            resolved_cutoffs[key] = float(val)

        if mol_names is not None:
            names_list = [str(n) for n in mol_names]
            if len(names_list) != len(self._ligand_entries):
                raise ValueError(
                    f"Length of mol_names ({len(names_list)}) must match number of ligands ({len(self._ligand_entries)})."
                )
        else:
            names_list = list(self.mol_names)

        per_mol_hits: list[dict[tuple[str, tuple], int]] = []
        per_mol_meta: list[tuple[str, str]] = []
        res_meta: dict[tuple, dict] = {}
        active_pairs: set[tuple[str, tuple]] = set()

        for mol_name, entry in zip(names_list, self._ligand_entries):
            mol_full = entry["ligand"]
            mol_no_h = entry["ligand_no_h"]
            h2f = entry["heavy_to_full_idx"]
            pocket_res = entry["pocket_residues"]
            smiles = Chem.MolToSmiles(mol_no_h)
            per_mol_meta.append((mol_name, smiles))
            for r in pocket_res:
                res_meta.setdefault(_residue_key(r), r)
            hits = self._count_prolig_interactions(
                mol_full=mol_full,
                mol_no_h=mol_no_h,
                heavy_to_full=h2f,
                pocket_residues=pocket_res,
                active_interactions=active,
                count_mode=count_mode,
                cutoffs=resolved_cutoffs,
            )
            per_mol_hits.append(hits)
            active_pairs.update(k for k, v in hits.items() if v > 0)

        if drop_empty:
            pairs_to_emit = active_pairs
            chain_pool = {rk for _, rk in active_pairs}
        else:
            pairs_to_emit = {
                (itype, rk)
                for itype in active
                for rk in res_meta
            }
            chain_pool = set(res_meta.keys())

        if show_chain == "auto":
            distinct_chains = {res_meta[rk]["chain"] for rk in chain_pool if res_meta[rk]["chain"]}
            include_chain = len(distinct_chains) > 1
        else:
            include_chain = bool(show_chain)

        order_idx = {name: i for i, name in enumerate(PROLIG_INTERACTION_ORDER)}
        sorted_pairs = sorted(
            pairs_to_emit,
            key=lambda item: (
                order_idx[item[0]],
                res_meta[item[1]]["residue_number"],
                res_meta[item[1]]["insertion_code"],
                res_meta[item[1]]["chain"],
                res_meta[item[1]]["residue_name"],
            ),
        )

        col_names: list[str] = []
        pair_to_col: dict[tuple[str, tuple], str] = {}
        seen_cols: set[str] = set()
        for itype, rk in sorted_pairs:
            col = f"{itype}_{_format_table_residue_tag(res_meta[rk], include_chain)}"
            pair_to_col[(itype, rk)] = col
            if col not in seen_cols:
                seen_cols.add(col)
                col_names.append(col)

        records = []
        for (mol_name, smiles), hits in zip(per_mol_meta, per_mol_hits):
            row: dict[str, Union[str, int]] = {"mol_name": mol_name, "smiles": smiles}
            for col in col_names:
                row[col] = 0
            for pair, val in hits.items():
                col = pair_to_col.get(pair)
                if col is not None and val > 0:
                    if count_mode == "binary":
                        row[col] = 1
                    else:
                        row[col] = int(row[col]) + int(val)
            records.append(row)

        df = pd.DataFrame(records, columns=["mol_name", "smiles", *col_names])
        for col in col_names:
            df[col] = df[col].astype(int)
        return df

    def plot_barcode(
        self,
        df: Optional[pd.DataFrame] = None,
        group_by: str = "interaction",
        show_counts: bool = False,
        legend: bool = True,
        color: Optional[dict] = None,
        xlabel: Optional[str] = None,
        dpi: int = 300,
        show: bool = True,
        savepath: Optional[str] = None,
        figsize: Optional[tuple[float, float]] = None,
        **to_df_kwargs,
    ):
        """
        Render a categorical protein–ligand interaction fingerprint barcode plot using Matplotlib.
        :param df: Optional[pd.DataFrame]
            Optional interaction count DataFrame produced by PLI.to_df(). When None (default), computed internally
            via ``self.to_df(**to_df_kwargs)``.
        :param group_by: str
            How Y-axis rows are grouped: 'interaction' (default, matches PLI.to_df() column order grouped by
            PROLIG_INTERACTION_ORDER with dividers between interaction types) or 'residue' (stacks interactions
            by residue in sequence order with dividers between residues).
        :param show_counts: bool
            If True, overlay integer contact counts centered inside occupied barcode cells. Defaults to False.
        :param legend: bool
            Whether to draw a compact interaction color legend to the right of the plot. Defaults to True.
        :param color: Optional[dict]
            Optional dictionary overriding default PROLIG_COLORS entries by interaction type (e.g.,
            {'HBDonor': '#1E40AF', 'Hydrophobic': '#2E8B57'}).
        :param xlabel: Optional[str]
            Optional X-axis title displayed below the molecule/pose labels.
        :param dpi: int
            Output resolution in dots per inch when saving raster figures. Defaults to 300.
        :param show: bool
            If True (default), calls plt.show() for inline display and returns None. If False, returns the
            Matplotlib Figure object.
        :param savepath: Optional[str]
            Optional output file path (.png, .svg, or any Matplotlib-supported vector/raster format).
        :param figsize: Optional[tuple[float, float]]
            Optional Matplotlib figure size (width, height) in inches. Auto-scaled from table dimensions when None.
        :param to_df_kwargs: dict
            Optional keyword arguments forwarded to PLI.to_df() when df is None (e.g., interactions='all',
            count='pairs'). When df is None, all loaded ligands are plotted automatically.
        :return: Optional[matplotlib.figure.Figure]
            Returns the Matplotlib Figure when show=False, otherwise None.
        """
        if group_by not in ("residue", "interaction"):
            raise ValueError(f"Invalid group_by={group_by!r}. Expected 'residue' or 'interaction'.")

        if figsize is not None:
            if (
                not isinstance(figsize, (tuple, list))
                or len(figsize) != 2
                or float(figsize[0]) <= 0
                or float(figsize[1]) <= 0
            ):
                raise ValueError(f"figsize must be a (width, height) tuple of positive numbers, got {figsize!r}.")
            plot_figsize = (float(figsize[0]), float(figsize[1]))
        else:
            plot_figsize = None

        if df is None:
            df = self.to_df(**to_df_kwargs)
        else:
            if to_df_kwargs:
                raise ValueError("Cannot specify both df and to_df keyword arguments.")
            if not isinstance(df, pd.DataFrame):
                raise TypeError("df must be a pandas DataFrame.")

        palette = _resolve_color_palette(PROLIG_COLORS, color)
        mol_labels = df["mol_name"].astype(str).tolist() if "mol_name" in df.columns else [str(i) for i in df.index]
        parsed, ytick_pos, ytick_labels, dividers = _prepare_barcode_rows(df, group_by)

        n_mols = len(df)
        n_rows = len(parsed)
        counts = (
            np.nan_to_num(df[[col for col, _, _ in parsed]].to_numpy(dtype=float).T, nan=0.0)
            if n_mols > 0 and n_rows > 0
            else np.empty((0, 0), dtype=float)
        )

        if n_mols == 0 or n_rows == 0 or not np.any(counts > 0):
            fig, ax = plt.subplots(figsize=plot_figsize or (5.0, 2.5), dpi=150)
            fig.patch.set_facecolor("white")
            ax.set_facecolor("white")
            ax.text(0.5, 0.5, "No interactions detected", ha="center", va="center", fontsize=11, color="#475569")
            ax.axis("off")
            return _finalize_matplotlib_figure(fig, savepath=savepath, dpi=dpi, show=show)

        present_set = {parsed[r][1] for r in range(n_rows) if np.any(counts[r] > 0)}
        present_interactions = [it for it in PROLIG_INTERACTION_ORDER if it in present_set]
        for _, itype, _ in parsed:
            if itype in present_set and itype not in present_interactions:
                present_interactions.append(itype)

        cat_id = {itype: idx + 1 for idx, itype in enumerate(present_interactions)}
        z_mat = np.zeros((n_rows, n_mols), dtype=int)
        for r, (_, itype, _) in enumerate(parsed):
            cid = cat_id.get(itype)
            if cid is not None:
                z_mat[r, counts[r] > 0] = cid

        cmap = ListedColormap(
            [(1.0, 1.0, 1.0)] + [palette.get(it, (0.45, 0.45, 0.45)) for it in present_interactions]
        )
        norm = BoundaryNorm(np.arange(-0.5, len(present_interactions) + 1.5, 1.0), cmap.N)

        if plot_figsize is None:
            height = max(2.8, min(18.0, 0.30 * n_rows + 1.2))
            base_w = max(2.4, min(14.0, 0.70 * n_mols + 1.2))
            legend_pad = 1.8 if legend else 0.4
            plot_figsize = (base_w + legend_pad, height)

        fig, ax = plt.subplots(figsize=plot_figsize, dpi=150)
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        ax.pcolormesh(
            np.arange(n_mols + 1),
            np.arange(n_rows + 1),
            z_mat,
            cmap=cmap,
            norm=norm,
            edgecolors="white",
            linewidth=0.8,
            shading="flat",
        )
        ax.set_xlim(0.0, float(n_mols))
        ax.set_ylim(float(n_rows), 0.0)

        for div_y in dividers:
            ax.axhline(y=div_y, color="#94A3B8", linewidth=0.9, zorder=3)

        ax.set_yticks(ytick_pos)
        ax.set_yticklabels(ytick_labels, fontsize=8.5, color="#1E293B")
        ax.tick_params(axis="y", length=0, pad=6)

        max_lbl_len = max((len(lbl) for lbl in mol_labels), default=0)
        if n_mols <= 5 and max_lbl_len <= 14:
            x_rot, x_ha, rot_mode = 0, "center", "default"
        elif n_mols <= 30:
            x_rot, x_ha, rot_mode = 45, "right", "anchor"
        else:
            x_rot, x_ha, rot_mode = 90, "center", "default"

        ax.set_xticks(np.arange(n_mols, dtype=float) + 0.5)
        tick_kw = {"fontsize": 8.5, "color": "#1E293B", "rotation": x_rot, "ha": x_ha}
        if rot_mode == "anchor":
            tick_kw["rotation_mode"] = "anchor"
        ax.set_xticklabels(mol_labels, **tick_kw)
        ax.tick_params(axis="x", length=3, color="#CBD5E1", pad=4)

        for spine in ax.spines.values():
            spine.set_edgecolor("#CBD5E1")
            spine.set_linewidth(0.8)

        if xlabel:
            ax.set_xlabel(str(xlabel), fontsize=9.5, color="#1E293B", labelpad=6)

        if show_counts:
            for r in range(n_rows):
                rgb = palette.get(parsed[r][1], (0.45, 0.45, 0.45))
                lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
                txt_color = "#1E293B" if lum > 0.62 else "#FFFFFF"
                for c in range(n_mols):
                    val = counts[r, c]
                    if val > 0:
                        ax.text(
                            c + 0.5,
                            r + 0.5,
                            f"{int(round(val))}",
                            ha="center",
                            va="center",
                            fontsize=7.8,
                            fontweight="normal",
                            color=txt_color,
                            zorder=4,
                        )

        if legend and present_interactions:
            handles = [
                patches.Patch(
                    facecolor=palette.get(it, (0.45, 0.45, 0.45)),
                    edgecolor="#94A3B8",
                    linewidth=0.5,
                    label=it,
                )
                for it in present_interactions
            ]
            ax.legend(
                handles=handles,
                bbox_to_anchor=(1.02, 1.0),
                loc="upper left",
                borderaxespad=0.0,
                frameon=True,
                facecolor="white",
                edgecolor="#CBD5E1",
                framealpha=0.95,
                fontsize=8.5,
                title="Interaction",
                title_fontsize=9.0,
            )

        fig.tight_layout()
        return _finalize_matplotlib_figure(fig, savepath=savepath, dpi=dpi, show=show)

    def _to_pharmacophore(self, **detect_kwargs) -> list:
        """
        Internal helper to convert detected 3D protein–ligand interactions into pharmacophore-toolkit's
        standard feature list.
        :param detect_kwargs: dict
            Optional keyword arguments forwarded to PLI._detect() (forces recomputation when supplied).
        :return: list
            List of [family: str, atom_indices: tuple, x: float, y: float, z: float] entries.
        """
        if detect_kwargs or self.interactions_df is None:
            self._detect(**detect_kwargs)

        lig_feats = self._classify_ligand_features()
        pos_atoms = set(lig_feats["pos_atoms"])

        features_list: list = []
        seen: set[tuple[str, frozenset[int]]] = set()
        aromatic_atoms: set[int] = set()
        hydrophobic_atoms: set[int] = set()

        for rec in self._interaction_records:
            itype = rec["interaction_type"]
            atom_indices = tuple(rec["ligand_atom_indices"])
            lx, ly, lz = rec["ligand_coords"]

            if itype == "Hydrophobic":
                hydrophobic_atoms.update(atom_indices)
                continue

            if itype == "HBond_Donor":
                family = "Donor"
            elif itype == "HBond_Acceptor":
                family = "Acceptor"
            elif itype == "Pi_Stacking":
                family = "Aromatic"
            elif itype == "Cation_Pi":
                family = "Aromatic" if len(atom_indices) > 1 else "PosIonizable"
            elif itype == "Ionic":
                family = "PosIonizable" if atom_indices[0] in pos_atoms else "NegIonizable"
            else:
                continue

            dedup_key = (family, frozenset(atom_indices))
            if dedup_key in seen:
                continue
            seen.add(dedup_key)

            if family == "Aromatic":
                aromatic_atoms.update(atom_indices)

            features_list.append(
                [family, atom_indices, round(float(lx), 4), round(float(ly), 4), round(float(lz), 4)]
            )

        remaining_hydrophobes = sorted(hydrophobic_atoms - aromatic_atoms)
        if remaining_hydrophobes:
            lig_conf = self.ligand_no_h.GetConformer()
            for patch in _connected_atom_patches(self.ligand_no_h, remaining_hydrophobes):
                coords = np.mean([lig_conf.GetAtomPosition(i) for i in patch], axis=0)
                features_list.append(
                    [
                        "Hydrophobe",
                        patch,
                        round(float(coords[0]), 4),
                        round(float(coords[1]), 4),
                        round(float(coords[2]), 4),
                    ]
                )

        self.features_list = features_list
        return features_list

    def output_pli(
        self,
        savepath: str,
        mol: Union[int, str] = 0,
        pocket_repr: Optional[str] = "lines",
        show_labels: bool = False,
        split_hbonds: bool = False,
        include_pharmacophore: bool = False,
        sphere_size: float = 0.5,
        line_width: float = 3.0,
        dash_gap: Optional[float] = None,
        dash_width: Optional[float] = 5.0,
        dash_radius: Optional[float] = 0.12,
        color: Optional[dict] = None,
        **detect_kwargs,
    ) -> str:
        """
        Write a drag-and-drop PyMOL (.pml) visualization script for the detected protein–ligand interactions.
        :param savepath: str
            Output file path ending in .pml.
        :param mol: Union[int, str]
            0-based index or molecule name specifying which loaded ligand to export. Defaults to ``0``.
        :param pocket_repr: Optional[str]
            PyMOL representation for interacting pocket residues ('lines' default, 'sticks', or None to omit).
        :param show_labels: bool
            Whether to show Ångström distance numbers on 3D dashed lines. Defaults to False.
        :param split_hbonds: bool
            If True, color HBond_Donor and HBond_Acceptor dashes separately; if False (default), use one unified
            HBond color.
        :param include_pharmacophore: bool
            If True, also emit 3D receptor-validated pharmacophore spheres from _to_pharmacophore().
        :param sphere_size: float
            PyMOL sphere_scale used when include_pharmacophore=True. Defaults to 0.5.
        :param line_width: float
            PyMOL line_width setting for pocket residues when pocket_repr='lines'. Defaults to 3.0.
        :param dash_gap: Optional[float]
            Optional PyMOL dash_gap override. Defaults to None (uses PyMOL's native Measurement Wizard setting).
        :param dash_width: Optional[float]
            PyMOL dash_width setting for interactive OpenGL dash thickness (None for PyMOL default). Defaults to 5.0.
        :param dash_radius: Optional[float]
            PyMOL dash_radius setting for ray-traced interaction dash thickness (None for PyMOL default). Defaults to 0.12.
        :param color: Optional[dict]
            Optional dictionary mapping interaction types (e.g., 'HBond', 'Pi_Stacking') to color names,
            hex strings, or RGB tuples.
        :param detect_kwargs: dict
            Optional keyword arguments forwarded to PLI._detect() (forces recomputation when supplied).
        :return: str
            Path to the written .pml file.
        """
        if not savepath or not isinstance(savepath, (str, os.PathLike)):
            raise ValueError("savepath must be a non-empty file path.")
        savepath_str = str(savepath).strip()
        if not savepath_str:
            raise ValueError("savepath must be a non-empty file path.")
        if pocket_repr not in ("lines", "sticks", None):
            raise ValueError(
                f"Invalid pocket_repr={pocket_repr!r}. Only 'lines', 'sticks', or None are supported."
            )

        self._select_ligand(mol)
        if detect_kwargs or self.interactions_df is None:
            self._detect(**detect_kwargs)

        palette = {k: list(v) for k, v in PML_COLORS.items()}
        if color is not None:
            if not isinstance(color, dict):
                raise ValueError("color must be a dictionary mapping interaction types to colors.")
            for key, val in color.items():
                if isinstance(val, (tuple, list)) and len(val) == 3:
                    palette[key] = [float(val[0]), float(val[1]), float(val[2])]
                else:
                    palette[key] = list(color_convert(val))

        parent_dir = os.path.dirname(savepath_str)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)

        dash_prefixes = {
            "HBond_Donor": "hbond",
            "HBond_Acceptor": "hbond",
            "Pi_Stacking": "pistack",
            "Cation_Pi": "catpi",
            "Ionic": "ionic",
            "Clash": "clash",
        }

        group_members: list[str] = []
        has_pocket = False

        with open(savepath_str, "w") as f:
            f.write("# Drag-and-drop PyMOL script generated by pharmacophore-toolkit (PLI)\n")
            f.write("set auto_zoom, off\n\n")

            for key, rgb in palette.items():
                f.write(f"set_color pli_{key.lower()}_color, [{rgb[0]:.4f}, {rgb[1]:.4f}, {rgb[2]:.4f}]\n")

            if include_pharmacophore:
                for fam, rgb in FEATURE_COLORS.items():
                    f.write(
                        f"set_color pli_pharm_{fam.lower()}_color, "
                        f"[{rgb[0]:.4f}, {rgb[1]:.4f}, {rgb[2]:.4f}]\n"
                    )
            f.write("recolor\n\n")

            unique_residues = sorted(
                {(rec["chain"], rec["residue_number"]) for rec in self._interaction_records}
            )
            if unique_residues:
                has_pocket = True
                clauses = [
                    f"(chain {ch} and resi {rnum})" if ch else f"(resi {rnum})"
                    for ch, rnum in unique_residues
                ]
                f.write(f"select pli_pocket_residues, polymer.protein and ({' or '.join(clauses)})\n")
                f.write("deselect\n")
                if pocket_repr is not None:
                    f.write(f"show {pocket_repr}, pli_pocket_residues\n")
                    if pocket_repr == "lines":
                        f.write(f"set line_width, {line_width}, pli_pocket_residues\n")
                    f.write("util.cnc pli_pocket_residues\n")
                f.write("\n")
                group_members.append("pli_pocket_residues")

            dash_idx = 0
            for rec in self._interaction_records:
                itype = rec["interaction_type"]
                if itype not in dash_prefixes:
                    continue
                dash_idx += 1
                prefix = dash_prefixes[itype]
                obj_name = f"{prefix}_{rec['residue_name']}{rec['residue_number']}_{dash_idx}"
                color_key = (
                    (itype if split_hbonds else "HBond")
                    if itype in ("HBond_Donor", "HBond_Acceptor")
                    else itype
                )
                lx, ly, lz = rec["ligand_coords"]
                px, py, pz = rec["protein_coords"]

                f.write(f"pseudoatom _pli_tmp_a, pos=[{lx}, {ly}, {lz}]\n")
                f.write(f"pseudoatom _pli_tmp_b, pos=[{px}, {py}, {pz}]\n")
                f.write(f"distance {obj_name}, _pli_tmp_a, _pli_tmp_b\n")
                f.write("delete _pli_tmp_a\n")
                f.write("delete _pli_tmp_b\n")
                f.write(f"set dash_color, pli_{color_key.lower()}_color, {obj_name}\n")
                if dash_gap is not None:
                    f.write(f"set dash_gap, {dash_gap}, {obj_name}\n")
                if dash_width is not None:
                    f.write(f"set dash_width, {dash_width}, {obj_name}\n")
                if dash_radius is not None:
                    f.write(f"set dash_radius, {dash_radius}, {obj_name}\n")
                if not show_labels:
                    f.write(f"hide labels, {obj_name}\n")
                f.write("\n")
                group_members.append(obj_name)

            if include_pharmacophore:
                feats = self._to_pharmacophore()
                coord_groups: dict[tuple[float, float, float], list[int]] = {}
                for idx, (_, _, x, y, z) in enumerate(feats):
                    key = (round(float(x), 3), round(float(y), 3), round(float(z), 3))
                    coord_groups.setdefault(key, []).append(idx)

                pharm_counts: dict[str, int] = {}
                split_specs = []
                for idx, (fam, _, x, y, z) in enumerate(feats):
                    pharm_counts[fam] = pharm_counts.get(fam, 0) + 1
                    obj_name = f"pli_pharm_{fam}_{pharm_counts[fam]}"
                    key = (round(float(x), 3), round(float(y), 3), round(float(z), 3))
                    group = coord_groups[key]
                    if len(group) == 1:
                        f.write(f"pseudoatom {obj_name}, pos=[{x}, {y}, {z}]\n")
                        f.write(f"show spheres, {obj_name}\n")
                        f.write(f"hide nonbonded, {obj_name}\n")
                        f.write(f"color pli_pharm_{fam.lower()}_color, {obj_name}\n")
                        f.write(f"set sphere_scale, {sphere_size}, {obj_name}\n")
                        f.write(f"set sphere_transparency, 0.25, {obj_name}\n\n")
                    else:
                        rank = group.index(idx)
                        k = len(group)
                        phi0 = 2.0 * math.pi * rank / k
                        phi1 = 2.0 * math.pi * (rank + 1) / k
                        rgb = tuple(FEATURE_COLORS.get(fam, (0.5, 0.5, 0.5)))
                        split_specs.append((obj_name, float(x), float(y), float(z), rgb, phi0, phi1))
                    group_members.append(obj_name)

                if split_specs:
                    f.write(
                        format_pymol_split_spheres_block(
                            split_specs,
                            sphere_scale=sphere_size,
                            transparency=0.25,
                        )
                    )

            if group_members:
                f.write(f"group pli_interactions, {' '.join(group_members)}\n")
            f.write("set auto_zoom, on\n")
            if has_pocket:
                f.write("zoom pli_pocket_residues, 6\n")

        print(f"Protein–ligand interaction script written to {savepath_str}.")
        return savepath_str

    def view(
        self,
        mol: Optional[Union[int, str]] = None,
        type: str = "jupyter",
        pocket_repr: str = "stick",
        protein_repr: Optional[str] = "cartoon",
        show_labels: bool = True,
        show_distances: bool = False,
        split_hbonds: bool = False,
        include_pharmacophore: bool = False,
        color: Optional[dict] = None,
        window: tuple[int, int] = (600, 500),
        prefix: str = "pli",
        **detect_kwargs,
    ):
        """
        Render an interactive 3D protein–ligand interaction view using py3Dmol (Jupyter) or 3Dmol.js (Marimo).
        :param mol: Optional[Union[int, str]]
            0-based index or molecule name specifying a single loaded ligand to display. When None (default),
            renders a single-ligand view if 1 ligand is loaded or an interactive dropdown selector across all loaded
            ligands if multiple ligands are present.
        :param type: str
            Notebook rendering backend: 'jupyter' (default) or 'marimo'.
        :param pocket_repr: str
            3D style for interacting pocket residues: 'stick' (default) or 'line'.
        :param protein_repr: Optional[str]
            3D style for the surrounding receptor chain: 'cartoon' (default), 'line', or None (pocket only).
        :param show_labels: bool
            Whether to display residue labels on contacting pocket residues. Defaults to True.
        :param show_distances: bool
            Whether to show Ångström distance badges at the midpoints of dashed interaction lines. Defaults to False.
        :param split_hbonds: bool
            If True, color HBond_Donor and HBond_Acceptor dashed cylinders separately; if False (default),
            use one unified HBond color.
        :param include_pharmacophore: bool
            If True, also render 3D receptor-validated pharmacophore spheres (with split wedges for coincident
            features). Defaults to False.
        :param color: Optional[dict]
            Optional dictionary mapping interaction types (e.g., 'HBond', 'Pi_Stacking') to color names,
            hex strings, or RGB tuples.
        :param window: tuple[int, int]
            Viewer (width, height) in pixels. Defaults to (600, 500).
        :param prefix: str
            Filename prefix for Marimo PNG snapshot exports. Defaults to 'pli'.
        :param detect_kwargs: dict
            Optional keyword arguments forwarded to PLI._detect() (forces recomputation when supplied).
        :return:
            Interactive py3Dmol viewer (or ipywidgets dropdown selector in Jupyter when multiple ligands are loaded)
            or Marimo iframe widget.
        """
        backend = str(type).lower()
        if backend not in ("jupyter", "marimo"):
            raise ValueError("Only 'jupyter' or 'marimo' accepted for type!")

        pocket_norm = {"stick": "stick", "sticks": "stick", "line": "line", "lines": "line"}.get(pocket_repr)
        if pocket_norm is None:
            raise ValueError(f"Invalid pocket_repr={pocket_repr!r}. Expected 'stick' or 'line'.")

        if protein_repr is None:
            prot_norm = None
        else:
            prot_norm = {"cartoon": "cartoon", "line": "line", "lines": "line"}.get(protein_repr)
            if prot_norm is None:
                raise ValueError(
                    f"Invalid protein_repr={protein_repr!r}. Expected 'cartoon', 'line', or None."
                )

        if not isinstance(window, (tuple, list)) or len(window) != 2 or int(window[0]) <= 0 or int(window[1]) <= 0:
            raise ValueError(f"window must be a (width, height) tuple of positive integers, got {window!r}.")
        win_size = (int(window[0]), int(window[1]))

        with open(self.protein_pdb, "r") as f:
            pdb_data = f.read()
        labels = _dropdown_labels(self.mol_names)

        def _build_payload(idx: int, include_pdb: bool = True) -> dict:
            self._select_ligand(idx)
            if detect_kwargs or self.interactions_df is None:
                self._detect(**detect_kwargs)
            return _prepare_view_payload(
                pli=self,
                pocket_repr=pocket_norm,
                protein_repr=prot_norm,
                show_labels=show_labels,
                show_distances=show_distances,
                split_hbonds=split_hbonds,
                include_pharmacophore=include_pharmacophore,
                color=color,
                pdb_data=pdb_data if include_pdb else None,
                include_pdb=include_pdb,
                name=labels[idx],
            )

        if mol is not None or len(self._ligand_entries) == 1:
            target_idx = 0 if mol is None else self._resolve_ligand_index(mol)
            if backend == "jupyter":
                return _render_pli_jupyter(_build_payload(target_idx, include_pdb=True), window=win_size)
            return _render_pli_marimo(
                pdb_data=pdb_data,
                lig_payloads=[_build_payload(target_idx, include_pdb=False)],
                protein_repr=prot_norm,
                pocket_repr=pocket_norm,
                window=win_size,
                prefix=prefix,
            )

        if backend == "jupyter":
            import ipywidgets as widgets

            dropdown = widgets.Dropdown(
                options=[(lbl, i) for i, lbl in enumerate(labels)],
                value=0,
                description="Select:",
                style={"description_width": "initial"},
            )

            def _render_selected(index: int):
                return _render_pli_jupyter(_build_payload(index, include_pdb=True), window=win_size)

            return widgets.interact(_render_selected, index=dropdown)

        prev_idx = self._active_idx
        try:
            lig_payloads = [
                _build_payload(idx, include_pdb=False)
                for idx in range(len(self._ligand_entries))
            ]
        finally:
            if prev_idx >= 0:
                self._select_ligand(prev_idx)
        return _render_pli_marimo(
            pdb_data=pdb_data,
            lig_payloads=lig_payloads,
            protein_repr=prot_norm,
            pocket_repr=pocket_norm,
            window=win_size,
            prefix=prefix,
        )

    def _prepare_2d_ligand(
        self,
        canvas_size: tuple[int, int] = (900, 900),
        padding: float = 0.28,
        black_and_white: bool = False,
    ) -> dict:
        """
        Compute a 2D depiction of the heavy-atom ligand on a padded Cairo canvas and map atom/ring pixel positions.
        :param canvas_size: tuple[int, int]
            Width and height of the 2D Cairo drawing canvas in pixels. Defaults to (900, 900).
        :param padding: float
            Fractional outer margin per side on the RDKit canvas (must satisfy 0.0 <= padding < 0.48).
            Defaults to 0.28.
        :param black_and_white: bool
            If True, override element atom colors with an all-black palette; if False (default), use RDKit's
            standard element colors.
        :return: dict
            Dictionary containing 'mol_2d', 'png_bytes', 'image', 'canvas_size', 'padding', 'atom_pixels',
            'ring_pixels', 'bond_segments', 'bond_pixel_length', and 'centroid'.
        """
        if not isinstance(canvas_size, (tuple, list)) or len(canvas_size) != 2:
            raise ValueError("canvas_size must be a (width, height) tuple of positive integers.")
        width, height = int(canvas_size[0]), int(canvas_size[1])
        if width <= 0 or height <= 0:
            raise ValueError(f"canvas_size dimensions must be positive, got {(width, height)}.")
        padding = float(padding)
        if not (0.0 <= padding < 0.48):
            raise ValueError(f"padding must satisfy 0.0 <= padding < 0.48, got {padding}.")

        # Copy self.ligand_no_h so the stored 3D conformers remain untouched
        mol_2d = Chem.Mol(self.ligand_no_h)
        rdDepictor.Compute2DCoords(mol_2d)

        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        opts = drawer.drawOptions()
        opts.padding = padding
        opts.bondLineWidth = 2.2
        if black_and_white:
            opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})

        drawer.DrawMolecule(mol_2d)
        drawer.FinishDrawing()
        png_bytes = drawer.GetDrawingText()
        image = np.array(Image.open(io.BytesIO(png_bytes)).convert("RGBA"))

        atom_pixels: dict[int, tuple[float, float]] = {}
        for idx in range(mol_2d.GetNumAtoms()):
            pt = drawer.GetDrawCoords(idx)
            atom_pixels[idx] = (float(pt.x), float(pt.y))

        ring_pixels: list[dict] = []
        for ring in mol_2d.GetRingInfo().AtomRings():
            coords = np.array([atom_pixels[i] for i in ring], dtype=float)
            cx, cy = coords.mean(axis=0)
            is_arom = all(mol_2d.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
            ring_pixels.append(
                {
                    "atom_indices": tuple(int(i) for i in ring),
                    "centroid": (float(cx), float(cy)),
                    "is_aromatic": bool(is_arom),
                }
            )

        bond_segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
        bond_lengths: list[float] = []
        for bond in mol_2d.GetBonds():
            p0 = atom_pixels[bond.GetBeginAtomIdx()]
            p1 = atom_pixels[bond.GetEndAtomIdx()]
            bond_segments.append((p0, p1))
            bond_lengths.append(float(np.hypot(p1[0] - p0[0], p1[1] - p0[1])))

        bond_pixel_length = float(np.median(bond_lengths)) if bond_lengths else 45.0
        if atom_pixels:
            pts = np.array(list(atom_pixels.values()), dtype=float)
            centroid = (float(pts[:, 0].mean()), float(pts[:, 1].mean()))
        else:
            centroid = (width / 2.0, height / 2.0)

        return {
            "mol_2d": mol_2d,
            "png_bytes": png_bytes,
            "image": image,
            "canvas_size": (width, height),
            "padding": padding,
            "atom_pixels": atom_pixels,
            "ring_pixels": ring_pixels,
            "bond_segments": bond_segments,
            "bond_pixel_length": bond_pixel_length,
            "centroid": centroid,
        }

    def draw(
        self,
        mol: Union[int, str] = 0,
        show_chain: Union[bool, str] = True,
        show_distances: bool = False,
        split_hbonds: bool = False,
        legend: bool = True,
        black_and_white: bool = False,
        color: Optional[dict] = None,
        dpi: int = 300,
        show: bool = True,
        savefig: Optional[str] = None,
        figsize: Optional[tuple[float, float]] = None,
        **detect_kwargs,
    ):
        """
        Render a 2D protein–ligand interaction diagram using RDKit and Matplotlib.
        :param mol: Union[int, str]
            0-based index or molecule name specifying which loaded ligand to draw. Defaults to ``0``.
        :param show_chain: Union[bool, str]
            Whether to append the PDB chain ID to residue labels (e.g., 'Phe97A'). Can be True (default), False,
            or 'auto' (shows chain IDs only when multiple chains contact the ligand).
        :param show_distances: bool
            Whether to draw Ångström distance pills at the midpoints of dashed connectors. Defaults to False.
        :param split_hbonds: bool
            If True, color HBond_Donor and HBond_Acceptor dashes separately; if False (default), use unified
            HBond color.
        :param legend: bool
            Whether to render a compact legend of present interaction types in the upper corner. Defaults to True.
        :param black_and_white: bool
            If True, render chemical atoms/labels in black & white while preserving interaction colors. Defaults to False.
        :param color: Optional[dict]
            Optional dictionary overriding default PLI_COLORS entries by interaction type (e.g.,
            {'Hydrophobic': '#4169E1'}).
        :param dpi: int
            Output resolution in dots per inch when saving raster figures. Defaults to 300.
        :param show: bool
            If True (default), calls plt.show() for inline display and returns None. If False, returns the
            Matplotlib Figure object.
        :param savefig: Optional[str]
            Optional output file path (.png, .svg, or any Matplotlib-supported vector/raster format).
        :param figsize: Optional[tuple[float, float]]
            Optional Matplotlib figure size (width, height) in inches. Auto-scaled from layout dimensions when None.
        :param detect_kwargs: dict
            Optional keyword arguments forwarded to PLI._detect() (forces re-detection when supplied). By
            default, draw() shows H-bonds, hydrophobic contacts, and aromatic (pi-pi / cation-pi) interactions
            with include_ionic=False; pass include_ionic=True to display ionic/salt bridges.
        :return: Optional[matplotlib.figure.Figure]
            Returns the Matplotlib Figure when show=False, otherwise None.
        """
        self._select_ligand(mol)
        draw_include_ionic = bool(detect_kwargs.get("include_ionic", False))
        if detect_kwargs:
            self._detect(**{"include_ionic": False, **detect_kwargs})
        elif self.interactions_df is None:
            self._detect(include_ionic=False)

        draw_records = [
            r
            for r in self._interaction_records
            if draw_include_ionic or r["interaction_type"] != "Ionic"
        ]

        palette = _resolve_color_palette(PLI_COLORS, color)

        if show_chain == "auto":
            chains = {r["chain"] for r in draw_records if r.get("chain")}
            include_chain = len(chains) > 1
        else:
            include_chain = bool(show_chain)

        scene = self._layout_poseview_scene(
            include_chain=include_chain,
            black_and_white=black_and_white,
            records=draw_records,
        )

        w, h = scene["canvas_size"]
        bond_px = scene["bond_pixel_length"]
        if figsize is None:
            figsize = (w / 120.0, h / 120.0)

        fig, ax = plt.subplots(figsize=figsize, dpi=150)
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        ax.imshow(scene["image"], extent=(0.0, float(w), float(h), 0.0), interpolation="bilinear")

        x_min, x_max = 0.0, float(w)
        y_min, y_max = 0.0, float(h)
        present_types: set[str] = set()

        for conn in scene["connectors"]:
            itype = conn["interaction_type"]
            color_key = (
                (itype if split_hbonds else "HBond")
                if itype in ("HBond_Donor", "HBond_Acceptor")
                else itype
            )
            line_color = palette.get(color_key, (0.2, 0.2, 0.2))
            present_types.add(color_key)

            p0 = np.array(conn["lig_point"], dtype=float)
            p1 = np.array(conn["prot_point"], dtype=float)

            if conn["kind"] == "aromatic":
                dot_r = 0.16 * bond_px
                if conn["lig_is_ring"]:
                    ax.add_patch(
                        patches.Circle(p0, radius=dot_r, facecolor=line_color, edgecolor="none", zorder=5)
                    )
                if conn["prot_is_ring"]:
                    ax.add_patch(
                        patches.Circle(p1, radius=dot_r, facecolor=line_color, edgecolor="none", zorder=5)
                    )
                s0 = p0
                if conn["prot_is_ring"]:
                    s1 = p1
                else:
                    s0, s1 = _clip_segment_for_glyphs(p0, p1, 0.0, conn["prot_pad_px"])
            else:
                s0, s1 = _clip_segment_for_glyphs(p0, p1, conn["lig_pad_px"], conn["prot_pad_px"])

            if itype in ("HBond_Donor", "HBond_Acceptor"):
                start, end = (s0, s1) if itype == "HBond_Donor" else (s1, s0)
                _draw_dashed_arrow(ax, start, end, line_color, bond_px)
            else:
                ax.add_line(
                    Line2D(
                        [s0[0], s1[0]],
                        [s0[1], s1[1]],
                        linestyle=(0, (4.2, 2.8)),
                        linewidth=1.9,
                        color=line_color,
                        solid_capstyle="butt",
                        zorder=4,
                    )
                )
            if show_distances and conn["distance"] is not None:
                mid = 0.5 * (s0 + s1)
                ax.text(
                    mid[0],
                    mid[1],
                    f"{conn['distance']:.1f} Å",
                    fontsize=7.5,
                    color=line_color,
                    ha="center",
                    va="center",
                    bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=line_color, lw=0.6, alpha=0.92),
                    zorder=7,
                )

        hydro_color = palette.get("Hydrophobic", PLI_COLORS["Hydrophobic"])
        for patch_info in scene["hydrophobic_patches"]:
            present_types.add("Hydrophobic")
            for arc_pts in patch_info["arcs"]:
                if len(arc_pts) < 2:
                    continue
                verts = np.asarray(arc_pts, dtype=float)
                ax.add_patch(
                    patches.PathPatch(
                        MplPath(verts),
                        facecolor="none",
                        edgecolor=hydro_color,
                        linewidth=1.85,
                        capstyle="round",
                        joinstyle="round",
                        zorder=3,
                    )
                )
                x_min = min(x_min, float(verts[:, 0].min()) - 12.0)
                x_max = max(x_max, float(verts[:, 0].max()) + 12.0)
                y_min = min(y_min, float(verts[:, 1].min()) - 12.0)
                y_max = max(y_max, float(verts[:, 1].max()) + 12.0)

            lx, ly = patch_info["label_pos"]
            ha = patch_info["label_ha"]
            labels = patch_info["labels"]
            label_text = "\n".join(labels)
            ax.text(
                lx,
                ly,
                label_text,
                color=hydro_color,
                fontsize=9.5,
                fontweight="normal",
                linespacing=1.25,
                ha=ha,
                va="center",
                zorder=6,
            )
            est_w, est_h = patch_info["label_box_size"]
            if ha == "left":
                bx0, bx1 = lx - 8.0, lx + est_w + 10.0
            elif ha == "right":
                bx0, bx1 = lx - est_w - 10.0, lx + 8.0
            else:
                bx0, bx1 = lx - 0.5 * est_w - 8.0, lx + 0.5 * est_w + 8.0
            x_min = min(x_min, bx0)
            x_max = max(x_max, bx1)
            y_min = min(y_min, ly - 0.5 * est_h - 10.0)
            y_max = max(y_max, ly + 0.5 * est_h + 10.0)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_max, y_min)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")

        if legend and present_types:
            _draw_pli_legend(ax, present_types, palette)

        return _finalize_matplotlib_figure(fig, savepath=savefig, dpi=dpi, show=show)

    def _layout_poseview_scene(
        self,
        include_chain: bool = True,
        black_and_white: bool = False,
        records: Optional[list[dict]] = None,
    ) -> dict:
        """
        Build a crossing-free 2D protein-ligand interaction multi-fragment scene (ligand + surrounding 2D residue
        fragments, explicit donor hydrogens, and hydrophobic contour splines).
        """
        if self.interactions_df is None:
            self._detect()

        active_records = self._interaction_records if records is None else records

        mol_2d = Chem.Mol(self.ligand_no_h)
        rdDepictor.Compute2DCoords(mol_2d)
        lig_conf = mol_2d.GetConformer()
        n_lig = mol_2d.GetNumAtoms()
        lig_pts = np.array(
            [[lig_conf.GetAtomPosition(i).x, lig_conf.GetAtomPosition(i).y] for i in range(n_lig)],
            dtype=float,
        )

        bond_lengths_2d = [
            float(np.linalg.norm(lig_pts[b.GetBeginAtomIdx()] - lig_pts[b.GetEndAtomIdx()]))
            for b in mol_2d.GetBonds()
        ]
        l_bond = float(np.median(bond_lengths_2d)) if bond_lengths_2d else 1.5
        lig_centroid = lig_pts.mean(axis=0) if len(lig_pts) else np.zeros(2, dtype=float)

        # Group directed contacts so multi-contact/bidentate residues share one 2D residue mini-structure
        grouped_records: dict[tuple, list[dict]] = {}
        for rec in active_records:
            itype = rec["interaction_type"]
            if itype == "Hydrophobic":
                continue
            is_arom = itype == "Pi_Stacking" or (
                itype == "Cation_Pi" and rec["protein_atom_name"] == "RING"
            )
            gkey = (
                rec["chain"],
                rec["residue_number"],
                rec["residue_name"],
                rec["is_backbone"],
                "arom" if is_arom else "polar",
            )
            grouped_records.setdefault(gkey, []).append(rec)

        frag_specs: list[dict] = []
        used_targets: list[np.ndarray] = []

        # Order around ligand consistently (aromatic/backbone/sidechain sorted deterministically)
        for gkey in sorted(grouped_records.keys(), key=lambda x: (x[0], x[1], x[4])):
            recs = grouped_records[gkey]
            recs_sorted = sorted(recs, key=lambda r: r["distance"])
            primary = recs_sorted[0]
            itype = primary["interaction_type"]
            lig_idxs = tuple(primary["ligand_atom_indices"])
            res_label = _format_poseview_label(
                primary["residue_name"],
                primary["residue_number"],
                primary["chain"],
                include_chain=include_chain,
            )
            tpl = _build_residue_2d_template(
                residue_name=primary["residue_name"],
                res_label=res_label,
                is_backbone=primary["is_backbone"],
                interaction_type=itype,
                protein_atom_name=primary["protein_atom_name"],
                target_bond_len=l_bond,
            )

            primary_xy = lig_pts[list(lig_idxs)].mean(axis=0)
            is_arom_group = gkey[4] == "arom"
            if len(recs_sorted) > 1 and not is_arom_group:
                anchor_xy = np.mean(
                    [lig_pts[list(r["ligand_atom_indices"])].mean(axis=0) for r in recs_sorted],
                    axis=0,
                )
            else:
                anchor_xy = primary_xy
            lig_is_ring = len(lig_idxs) > 1
            lig_is_donor = (itype == "HBond_Donor") and (len(lig_idxs) == 1)

            if lig_is_ring:
                theta_0 = float(np.arctan2(anchor_xy[1] - lig_centroid[1], anchor_xy[0] - lig_centroid[0]))
            else:
                a_idx = lig_idxs[0]
                nbrs = [nb.GetIdx() for nb in mol_2d.GetAtomWithIdx(a_idx).GetNeighbors()]
                if nbrs:
                    out_vec = np.zeros(2, dtype=float)
                    for nb_idx in nbrs:
                        diff = lig_pts[a_idx] - lig_pts[nb_idx]
                        nrm = float(np.linalg.norm(diff))
                        if nrm > 1e-8:
                            out_vec += diff / nrm
                    if np.linalg.norm(out_vec) > 1e-6:
                        theta_0 = float(np.arctan2(out_vec[1], out_vec[0]))
                    else:
                        theta_0 = float(np.arctan2(anchor_xy[1] - lig_centroid[1], anchor_xy[0] - lig_centroid[0]))
                else:
                    theta_0 = float(np.arctan2(anchor_xy[1] - lig_centroid[1], anchor_xy[0] - lig_centroid[0]))

            if is_arom_group and lig_is_ring:
                d_0 = 3.5 * l_bond
            elif lig_is_donor:
                d_0 = 3.3 * l_bond
            elif itype == "Cation_Pi":
                d_0 = 2.9 * l_bond
            else:
                d_0 = 2.5 * l_bond

            # Evaluate 72 initial outward directions for max clearance from ligand and earlier residues
            thetas = np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False)
            best_score = -1e18
            best_theta = theta_0
            for th in thetas:
                u = np.array([np.cos(th), np.sin(th)], dtype=float)
                probe1 = anchor_xy + 1.8 * l_bond * u
                probe2 = anchor_xy + d_0 * u
                other_mask = [i for i in range(n_lig) if i not in lig_idxs]
                ref_pts = lig_pts[other_mask] if other_mask else lig_pts
                min_d1 = float(np.linalg.norm(ref_pts - probe1, axis=1).min())
                min_d2 = float(np.linalg.norm(lig_pts - probe2, axis=1).min())
                score = min(min_d1, 1.5 * l_bond) + min_d2
                if used_targets:
                    sep = min(float(np.linalg.norm(probe2 - prev_q)) for prev_q in used_targets)
                    score += min(sep, 3.0 * l_bond) * 0.7
                if not lig_is_ring:
                    score += 0.9 * l_bond * float(np.cos(th - theta_0))
                if score > best_score:
                    best_score = score
                    best_theta = float(th)

            q_init = anchor_xy + d_0 * np.array([np.cos(best_theta), np.sin(best_theta)], dtype=float)
            used_targets.append(q_init)

            frag_specs.append(
                {
                    "records": recs_sorted,
                    "primary": primary,
                    "tpl": tpl,
                    "lig_idxs": lig_idxs,
                    "primary_xy": primary_xy,
                    "anchor_xy": anchor_xy,
                    "theta_init": best_theta,
                    "theta_pref": best_theta if lig_is_ring else theta_0,
                    "theta": best_theta,
                    "d_0": d_0,
                    "d": d_0,
                    "lig_is_ring": lig_is_ring,
                    "lig_is_donor": lig_is_donor,
                    "is_arom_group": is_arom_group,
                }
            )

        if frag_specs:
            _relax_2d_fragments(frag_specs, lig_pts, mol_2d, l_bond)

        # Orient final fragment coordinates and compute outgoing vectors
        for spec in frag_specs:
            pts, q_pos = _orient_fragment_coords(
                spec["tpl"],
                spec["anchor_xy"],
                spec["theta"],
                spec["d"],
                lig_pts,
                l_bond,
            )
            if (
                len(spec["records"]) > 1
                and spec["tpl"]["contact_idx"] is not None
                and spec["tpl"]["anchor_heavy_idx"] is not None
            ):
                c_idx = spec["tpl"]["contact_idx"]
                h_idx = spec["tpl"]["anchor_heavy_idx"]
                if c_idx != h_idx and spec["tpl"]["mol"].GetAtomWithIdx(c_idx).GetAtomicNum() == 1:
                    v_h = spec["primary_xy"] - pts[h_idx]
                    vn = float(np.linalg.norm(v_h))
                    if vn > 1e-8:
                        pts = pts.copy()
                        pts[c_idx] = pts[h_idx] + l_bond * (v_h / vn)
                        q_pos = pts[c_idx]
            spec["coords_2d"] = pts
            spec["prot_anchor_2d"] = q_pos

        # Sprout collinear explicit polar donor H atoms on mol_2d for ligand HBond_Donor contacts
        rw_lig = Chem.RWMol(mol_2d)
        sprout_positions: list[tuple[int, np.ndarray]] = []
        remaining_hs = {i: mol_2d.GetAtomWithIdx(i).GetTotalNumHs() for i in range(n_lig)}
        atom_sprouted_hs: dict[int, list[int]] = {}
        record_lig_atom_override: dict[int, int] = {}

        for spec_idx, spec in enumerate(frag_specs):
            for rec in spec["records"]:
                if rec["interaction_type"] != "HBond_Donor" or len(rec["ligand_atom_indices"]) != 1:
                    continue
                lig_a = rec["ligand_atom_indices"][0]
                if remaining_hs.get(lig_a, 0) > 0:
                    atom_obj = rw_lig.GetAtomWithIdx(lig_a)
                    if atom_obj.GetNumExplicitHs() > 0:
                        atom_obj.SetNumExplicitHs(atom_obj.GetNumExplicitHs() - 1)
                    h_idx = rw_lig.AddAtom(Chem.Atom(1))
                    rw_lig.AddBond(lig_a, h_idx, Chem.BondType.SINGLE)
                    remaining_hs[lig_a] -= 1

                    out_vec = spec["prot_anchor_2d"] - lig_pts[lig_a]
                    nrm = float(np.linalg.norm(out_vec))
                    u = out_vec / nrm if nrm > 1e-8 else np.array([1.0, 0.0])
                    h_pos = lig_pts[lig_a] + l_bond * u
                    sprout_positions.append((h_idx, h_pos))
                    atom_sprouted_hs.setdefault(lig_a, []).append(h_idx)
                    record_lig_atom_override[id(rec)] = h_idx
                elif atom_sprouted_hs.get(lig_a):
                    record_lig_atom_override[id(rec)] = atom_sprouted_hs[lig_a][0]

        mol_2d_final = rw_lig.GetMol()
        if sprout_positions:
            Chem.SanitizeMol(mol_2d_final)
            conf_final = mol_2d_final.GetConformer()
            for h_idx, h_pos in sprout_positions:
                conf_final.SetAtomPosition(h_idx, Point3D(float(h_pos[0]), float(h_pos[1]), 0.0))

        # Write oriented coordinates back into each fragment Mol and combine into one unified scene
        combo_mol = Chem.Mol(mol_2d_final)
        frag_offsets: list[int] = []
        for spec in frag_specs:
            fmol = Chem.Mol(spec["tpl"]["mol"])
            fconf = fmol.GetConformer()
            fpts = spec["coords_2d"]
            for i in range(fmol.GetNumAtoms()):
                fconf.SetAtomPosition(i, Point3D(float(fpts[i, 0]), float(fpts[i, 1]), 0.0))
            offset = combo_mol.GetNumAtoms()
            frag_offsets.append(offset)
            combo_mol = Chem.CombineMols(combo_mol, fmol)

        # Compute bounding box expanded around hydrophobic ligand atoms for ample margin
        combo_conf = combo_mol.GetConformer()
        all_xy = np.array(
            [
                [combo_conf.GetAtomPosition(i).x, combo_conf.GetAtomPosition(i).y]
                for i in range(combo_mol.GetNumAtoms())
            ],
            dtype=float,
        )
        hydro_lig_idxs = sorted(
            {
                idx
                for r in active_records
                if r["interaction_type"] == "Hydrophobic"
                for idx in r["ligand_atom_indices"]
            }
        )
        if hydro_lig_idxs:
            hpts = lig_pts[hydro_lig_idxs]
            min_x = min(float(all_xy[:, 0].min()), float(hpts[:, 0].min()) - 2.3 * l_bond)
            max_x = max(float(all_xy[:, 0].max()), float(hpts[:, 0].max()) + 2.3 * l_bond)
            min_y = min(float(all_xy[:, 1].min()), float(hpts[:, 1].min()) - 2.3 * l_bond)
            max_y = max(float(all_xy[:, 1].max()), float(hpts[:, 1].max()) + 2.3 * l_bond)
        else:
            min_x, max_x = float(all_xy[:, 0].min()), float(all_xy[:, 0].max())
            min_y, max_y = float(all_xy[:, 1].min()), float(all_xy[:, 1].max())

        px_per_unit = 54.0 / l_bond
        padding = 0.16
        width = int(np.clip((max_x - min_x) * px_per_unit / (1.0 - 2.0 * padding), 900, 1650))
        height = int(np.clip((max_y - min_y) * px_per_unit / (1.0 - 2.0 * padding), 680, 1400))

        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        opts = drawer.drawOptions()
        opts.padding = padding
        opts.bondLineWidth = 2.1
        if black_and_white:
            opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})
        else:
            opts.updateAtomPalette({0: (0.0, 0.0, 0.0), 1: (0.0, 0.0, 0.0)})

        drawer.DrawMolecule(combo_mol)
        drawer.FinishDrawing()
        png_bytes = drawer.GetDrawingText()
        image = np.array(Image.open(io.BytesIO(png_bytes)).convert("RGBA"))

        atom_pixels: dict[int, tuple[float, float]] = {}
        for idx in range(combo_mol.GetNumAtoms()):
            pt = drawer.GetDrawCoords(idx)
            atom_pixels[idx] = (float(pt.x), float(pt.y))

        lig_bond_px = [
            float(
                np.hypot(
                    atom_pixels[b.GetBeginAtomIdx()][0] - atom_pixels[b.GetEndAtomIdx()][0],
                    atom_pixels[b.GetBeginAtomIdx()][1] - atom_pixels[b.GetEndAtomIdx()][1],
                )
            )
            for b in mol_2d_final.GetBonds()
        ]
        bond_px = float(np.median(lig_bond_px)) if lig_bond_px else 50.0

        # Build pixel-coordinate connectors
        connectors: list[dict] = []
        polar_lig_heavy_atoms: set[int] = set()
        for spec, offset in zip(frag_specs, frag_offsets):
            tpl = spec["tpl"]
            for rec in spec["records"]:
                itype = rec["interaction_type"]
                l_idxs = tuple(rec["ligand_atom_indices"])
                if len(l_idxs) == 1:
                    polar_lig_heavy_atoms.add(l_idxs[0])

                is_arom = itype == "Pi_Stacking" or (
                    itype == "Cation_Pi" and rec["protein_atom_name"] == "RING"
                )
                lig_is_ring = len(l_idxs) > 1
                if lig_is_ring:
                    l_xy = tuple(np.mean([atom_pixels[i] for i in l_idxs], axis=0))
                    lig_pad_px = 0.0
                else:
                    l_atom_idx = record_lig_atom_override.get(id(rec), l_idxs[0])
                    l_xy = atom_pixels[l_atom_idx]
                    lig_pad_px = _atom_glyph_radius_px(mol_2d_final.GetAtomWithIdx(l_atom_idx), bond_px)

                prot_ring_idxs = tpl["ring_indices"]
                prot_is_ring = is_arom and prot_ring_idxs is not None
                if prot_is_ring:
                    p_xy = tuple(np.mean([atom_pixels[offset + i] for i in prot_ring_idxs], axis=0))
                    prot_pad_px = 0.0
                else:
                    p_local_idx = _resolve_template_contact_idx(tpl, rec)
                    p_global_idx = offset + p_local_idx
                    p_xy = atom_pixels[p_global_idx]
                    prot_pad_px = _atom_glyph_radius_px(combo_mol.GetAtomWithIdx(p_global_idx), bond_px)

                connectors.append(
                    {
                        "interaction_type": itype,
                        "kind": "aromatic" if is_arom else "polar",
                        "lig_point": l_xy,
                        "prot_point": p_xy,
                        "lig_is_ring": lig_is_ring,
                        "prot_is_ring": prot_is_ring,
                        "lig_pad_px": lig_pad_px,
                        "prot_pad_px": prot_pad_px,
                        "distance": rec["distance"],
                    }
                )

        # Build pixel-space hydrophobic patches, contour splines, and stacked unboxed labels
        frag_atom_pixels = [atom_pixels[i] for i in range(mol_2d_final.GetNumAtoms(), combo_mol.GetNumAtoms())]
        hydrophobic_patches = _build_hydrophobic_splines_and_labels(
            records=active_records,
            mol_2d=mol_2d_final,
            n_heavy_lig=n_lig,
            atom_pixels=atom_pixels,
            frag_atom_pixels=frag_atom_pixels,
            connectors=connectors,
            polar_lig_heavy_atoms=polar_lig_heavy_atoms,
            bond_px=bond_px,
            canvas_size=(width, height),
            include_chain=include_chain,
        )

        return {
            "combo_mol": combo_mol,
            "image": image,
            "canvas_size": (width, height),
            "atom_pixels": atom_pixels,
            "bond_pixel_length": bond_px,
            "connectors": connectors,
            "hydrophobic_patches": hydrophobic_patches,
        }

    def _classify_ligand_features(
        self,
        lig_full: Optional[Chem.Mol] = None,
        lig_no_h: Optional[Chem.Mol] = None,
        heavy_to_full: Optional[dict[int, int]] = None,
    ) -> dict:
        """
        Identify donors, acceptors, halogen donors, metal centers/acceptors, ionizable groups, aromatic
        rings, and hydrophobic atoms on a ligand.
        """
        lig_full = self.ligand if lig_full is None else lig_full
        lig_no_h = self.ligand_no_h if lig_no_h is None else lig_no_h
        heavy_to_full = self.heavy_to_full_idx if heavy_to_full is None else heavy_to_full

        donors: list[tuple[int, int, tuple[int, ...]]] = []
        acceptors: list[tuple[int, int]] = []
        halogen_donors: list[tuple[int, int]] = []
        metal_atoms: list[int] = []
        metal_acceptors: list[int] = []
        hydrophobes: set[int] = set()

        for heavy_idx in range(lig_no_h.GetNumAtoms()):
            full_idx = heavy_to_full[heavy_idx]
            full_atom = lig_full.GetAtomWithIdx(full_idx)
            heavy_atom = lig_no_h.GetAtomWithIdx(heavy_idx)
            anum = heavy_atom.GetAtomicNum()
            chg = heavy_atom.GetFormalCharge()

            if anum in (7, 8, 16):
                explicit_hs = tuple(
                    nbr.GetIdx() for nbr in full_atom.GetNeighbors() if nbr.GetAtomicNum() == 1
                )
                if explicit_hs or heavy_atom.GetTotalNumHs() > 0:
                    donors.append((heavy_idx, full_idx, explicit_hs))

                if anum == 8 and chg <= 0:
                    acceptors.append((heavy_idx, full_idx))
                elif anum == 7 and is_hbond_acceptor_nitrogen(full_atom):
                    acceptors.append((heavy_idx, full_idx))
                elif anum == 16 and chg <= 0 and heavy_atom.GetDegree() <= 2:
                    acceptors.append((heavy_idx, full_idx))

            if anum in METAL_ATOMIC_NUMS:
                metal_atoms.append(heavy_idx)
            elif chg <= 0 and (anum in (7, 8, 16) or chg < 0):
                metal_acceptors.append(heavy_idx)

            if anum in (17, 35, 53):
                c_nbrs = [nb.GetIdx() for nb in heavy_atom.GetNeighbors() if nb.GetAtomicNum() == 6]
                if c_nbrs:
                    halogen_donors.append((heavy_idx, c_nbrs[0]))

            is_nonpolar_c = (
                anum == 6
                and chg == 0
                and not any(nb.GetAtomicNum() in (7, 8, 15, 16) for nb in heavy_atom.GetNeighbors())
            )
            if is_nonpolar_c or anum in (9, 17, 35, 53):
                hydrophobes.add(heavy_idx)

        pos_atoms: set[int] = set()
        for q in (_POS_CHARGED_N_QUERY, _POS_AMINE_QUERY):
            if q is not None:
                for match in lig_no_h.GetSubstructMatches(q):
                    pos_atoms.add(match[0])

        neg_atoms: set[int] = {
            atom.GetIdx()
            for atom in lig_no_h.GetAtoms()
            if atom.GetFormalCharge() < 0
            and not any(nb.GetFormalCharge() > 0 for nb in atom.GetNeighbors())
        }
        if _NEG_ACID_O_QUERY is not None:
            for match in lig_no_h.GetSubstructMatches(_NEG_ACID_O_QUERY):
                neg_atoms.add(match[0])

        lig_conf = lig_no_h.GetConformer()
        arom_rings: list[tuple[tuple[int, ...], np.ndarray, np.ndarray]] = []
        for ring in lig_no_h.GetRingInfo().AtomRings():
            if len(ring) in (5, 6) and all(
                lig_no_h.GetAtomWithIdx(i).GetIsAromatic() for i in ring
            ):
                coords = np.array([lig_conf.GetAtomPosition(i) for i in ring], dtype=float)
                centroid = coords.mean(axis=0)
                normal = _ring_normal(coords)
                arom_rings.append((tuple(ring), centroid, normal))

        return {
            "donors": donors,
            "acceptors": acceptors,
            "halogen_donors": halogen_donors,
            "metal_atoms": metal_atoms,
            "metal_acceptors": metal_acceptors,
            "pos_atoms": sorted(pos_atoms),
            "neg_atoms": sorted(neg_atoms),
            "arom_rings": arom_rings,
            "hydrophobes": hydrophobes,
        }

    def _detect_hbonds(
        self,
        lig_feats: dict,
        hbond_dist: float,
        hbond_angle: float,
    ) -> tuple[list[dict], set[tuple[tuple, int]]]:
        """
        Detect hydrogen bonds between pocket residues and ligand donors/acceptors.
        """
        prot_conf = self.protein.GetConformer()
        lig_conf_full = self.ligand.GetConformer()

        prot_donors = []
        prot_acceptors = []
        for res in self.pocket_residues:
            res_key = _residue_key(res)
            for p_idx in res["heavy_atom_indices"]:
                p_atom = self.protein.GetAtomWithIdx(p_idx)
                if p_atom.GetAtomicNum() not in (7, 8, 16):
                    continue
                info = p_atom.GetPDBResidueInfo()
                atom_name = info.GetName().strip().upper()
                can_don, can_acc, is_bb = _protein_hbond_role(
                    res["residue_name"], res["pdb_residue_name"], atom_name, p_atom
                )
                p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                if can_don:
                    prot_donors.append((res, res_key, p_idx, atom_name, is_bb, p_atom, p_pos))
                if can_acc:
                    prot_acceptors.append((res, res_key, p_idx, atom_name, is_bb, p_atom, p_pos))

        best: dict[tuple[tuple, str, int], dict] = {}

        # Ligand Donor -> Protein Acceptor
        for lig_h_idx, lig_f_idx, _ in lig_feats["donors"]:
            d_atom = self.ligand.GetAtomWithIdx(lig_f_idx)
            d_pos = np.array(lig_conf_full.GetAtomPosition(lig_f_idx), dtype=float)
            for res, res_key, _, atom_name, is_bb, a_atom, a_pos in prot_acceptors:
                dist = float(np.linalg.norm(d_pos - a_pos))
                if not (2.0 <= dist <= hbond_dist):
                    continue
                ok, ang = _eval_hbond_geometry(
                    d_atom,
                    d_pos,
                    a_pos,
                    lig_conf_full,
                    hbond_angle,
                    DEFAULT_HBOND_HEAVY_ANGLE,
                    acceptor_atom=a_atom,
                    acceptor_conf=prot_conf,
                )
                if not ok:
                    continue
                key = (res_key, "HBond_Donor", lig_h_idx)
                if key not in best or dist < best[key]["distance"]:
                    best[key] = _make_row(
                        res, res_key, is_bb, "HBond_Donor", (lig_h_idx,),
                        atom_name, dist, ang, d_pos, a_pos,
                    )

        # Protein Donor -> Ligand Acceptor
        for res, res_key, _, atom_name, is_bb, d_atom, d_pos in prot_donors:
            for lig_h_idx, lig_f_idx in lig_feats["acceptors"]:
                a_atom = self.ligand.GetAtomWithIdx(lig_f_idx)
                a_pos = np.array(lig_conf_full.GetAtomPosition(lig_f_idx), dtype=float)
                dist = float(np.linalg.norm(d_pos - a_pos))
                if not (2.0 <= dist <= hbond_dist):
                    continue
                ok, ang = _eval_hbond_geometry(
                    d_atom,
                    d_pos,
                    a_pos,
                    prot_conf,
                    hbond_angle,
                    DEFAULT_HBOND_HEAVY_ANGLE,
                    acceptor_atom=a_atom,
                    acceptor_conf=lig_conf_full,
                )
                if not ok:
                    continue
                key = (res_key, "HBond_Acceptor", lig_h_idx)
                if key not in best or dist < best[key]["distance"]:
                    best[key] = _make_row(
                        res, res_key, is_bb, "HBond_Acceptor", (lig_h_idx,),
                        atom_name, dist, ang, a_pos, d_pos,
                    )

        rows = list(best.values())
        hbond_pairs = {(r["_res_key"], r["ligand_atom_indices"][0]) for r in rows}
        return rows, hbond_pairs

    def _detect_ionic(
        self,
        lig_feats: dict,
        ionic_dist: float,
        hbond_pairs: set[tuple[tuple, int]],
    ) -> list[dict]:
        """
        Detect salt bridges / ionic interactions between oppositely charged ligand and protein atoms.
        """
        prot_conf = self.protein.GetConformer()
        lig_conf = self.ligand_no_h.GetConformer()
        best: dict[tuple[tuple, int], dict] = {}

        for res in self.pocket_residues:
            res_key = _residue_key(res)
            pos_names, neg_names = _protein_pos_neg_atom_names(res)
            if not pos_names and not neg_names:
                continue

            for p_idx in res["heavy_atom_indices"]:
                p_atom = self.protein.GetAtomWithIdx(p_idx)
                atom_name = p_atom.GetPDBResidueInfo().GetName().strip().upper()
                p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)

                target_lig_atoms: list[int] = []
                if atom_name in pos_names:
                    target_lig_atoms = lig_feats["neg_atoms"]
                elif atom_name in neg_names:
                    target_lig_atoms = lig_feats["pos_atoms"]

                for l_idx in target_lig_atoms:
                    if (res_key, l_idx) in hbond_pairs:
                        continue
                    l_pos = np.array(lig_conf.GetAtomPosition(l_idx), dtype=float)
                    dist = float(np.linalg.norm(p_pos - l_pos))
                    if dist <= ionic_dist:
                        key = (res_key, l_idx)
                        if key not in best or dist < best[key]["distance"]:
                            best[key] = _make_row(
                                res, res_key, False, "Ionic", (l_idx,),
                                atom_name, dist, None, l_pos, p_pos,
                            )

        return list(best.values())

    def _detect_aromatic(
        self,
        lig_feats: dict,
        pi_stacking_dist: float,
        cation_pi_dist: float,
        pi_face_dist: float = DEFAULT_PI_FACE_DIST,
        pi_face_angle: float = DEFAULT_PI_FACE_ANGLE,
        pi_edge_angle: float = DEFAULT_PI_EDGE_ANGLE,
        cation_pi_angle: float = DEFAULT_CATION_PI_ANGLE,
    ) -> list[dict]:
        """
        Detect aromatic pi-pi stacking and cation-pi interactions.
        """
        prot_conf = self.protein.GetConformer()
        lig_conf = self.ligand_no_h.GetConformer()

        best_pi: dict[tuple[tuple, tuple[int, ...]], dict] = {}
        best_cat: dict[tuple[tuple, tuple[int, ...]], dict] = {}

        for res in self.pocket_residues:
            res_key = _residue_key(res)
            canon = res["residue_name"]

            # Protein aromatic rings -> Ligand aromatic rings (Pi_Stacking) and Ligand cations (Cation_Pi)
            if canon in PROTEIN_AROMATIC_RESIDUES:
                prot_rings = self._extract_protein_aromatic_rings(res)
                for _, p_cent, p_norm in prot_rings:
                    for l_ring, l_cent, l_norm in lig_feats["arom_rings"]:
                        dist = float(np.linalg.norm(p_cent - l_cent))
                        if dist > pi_stacking_dist:
                            continue
                        theta = float(
                            np.degrees(np.arccos(np.clip(abs(float(np.dot(p_norm, l_norm))), 0.0, 1.0)))
                        )
                        is_parallel = dist <= pi_face_dist and theta <= pi_face_angle
                        is_t_shaped = dist <= pi_stacking_dist and pi_edge_angle <= theta <= 90.0
                        if is_parallel or is_t_shaped:
                            key = (res_key, l_ring)
                            if key not in best_pi or dist < best_pi[key]["distance"]:
                                best_pi[key] = _make_row(
                                    res, res_key, False, "Pi_Stacking", l_ring,
                                    "RING", dist, theta, l_cent, p_cent,
                                )

                    for l_pos_idx in lig_feats["pos_atoms"]:
                        l_pos = np.array(lig_conf.GetAtomPosition(l_pos_idx), dtype=float)
                        dist = float(np.linalg.norm(p_cent - l_pos))
                        if dist <= cation_pi_dist and _cation_pi_angle(l_pos, p_cent, p_norm) <= cation_pi_angle:
                            key = (res_key, (l_pos_idx,))
                            if key not in best_cat or dist < best_cat[key]["distance"]:
                                best_cat[key] = _make_row(
                                    res, res_key, False, "Cation_Pi", (l_pos_idx,),
                                    "RING", dist, None, l_pos, p_cent,
                                )

            # Protein cations (LYS NZ / ARG CZ / HIP ND1,NE2) -> Ligand aromatic rings (Cation_Pi)
            cation_names = _protein_pi_cation_names(res)
            if cation_names and lig_feats["arom_rings"]:
                for p_idx in res["heavy_atom_indices"]:
                    atom_name = self.protein.GetAtomWithIdx(p_idx).GetPDBResidueInfo().GetName().strip().upper()
                    if atom_name not in cation_names:
                        continue
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    for l_ring, l_cent, l_norm in lig_feats["arom_rings"]:
                        dist = float(np.linalg.norm(p_pos - l_cent))
                        if dist <= cation_pi_dist and _cation_pi_angle(p_pos, l_cent, l_norm) <= cation_pi_angle:
                            key = (res_key, l_ring)
                            if key not in best_cat or dist < best_cat[key]["distance"]:
                                best_cat[key] = _make_row(
                                    res, res_key, False, "Cation_Pi", l_ring,
                                    atom_name, dist, None, l_cent, p_pos,
                                )

        return list(best_pi.values()) + list(best_cat.values())

    def _extract_protein_aromatic_rings(
        self, res_entry: dict
    ) -> list[tuple[tuple[int, ...], np.ndarray, np.ndarray]]:
        """
        Extract (atom_indices, centroid, unit_normal) for aromatic sidechain rings of a pocket residue.
        """
        prot_conf = self.protein.GetConformer()
        res_set = set(res_entry["heavy_atom_indices"])
        seen: set[frozenset[int]] = set()
        rings: list[tuple[tuple[int, ...], np.ndarray, np.ndarray]] = []

        for ring in self.protein.GetRingInfo().AtomRings():
            if len(ring) in (5, 6) and set(ring).issubset(res_set):
                fset = frozenset(ring)
                seen.add(fset)
                coords = np.array([prot_conf.GetAtomPosition(i) for i in ring], dtype=float)
                rings.append((tuple(ring), coords.mean(axis=0), _ring_normal(coords)))

        # Fallback by standard PDB atom names when sanitization/bond order skipped a ring
        name_to_idx = {
            self.protein.GetAtomWithIdx(i).GetPDBResidueInfo().GetName().strip().upper(): i
            for i in res_entry["heavy_atom_indices"]
        }
        for ring_names in PROTEIN_AROMATIC_RING_ATOMS.get(res_entry["residue_name"], ()):
            if all(n in name_to_idx for n in ring_names):
                ring = tuple(name_to_idx[n] for n in ring_names)
                fset = frozenset(ring)
                if fset not in seen:
                    seen.add(fset)
                    coords = np.array([prot_conf.GetAtomPosition(i) for i in ring], dtype=float)
                    rings.append((ring, coords.mean(axis=0), _ring_normal(coords)))

        return rings

    def _detect_hydrophobic(
        self,
        lig_feats: dict,
        hydrophobic_dist: float,
        specific_contact_residues: set[tuple],
    ) -> list[dict]:
        """
        Detect hydrophobic contacts grouped into one row per contacting hydrophobic residue.
        """
        if not lig_feats["hydrophobes"]:
            return []

        prot_conf = self.protein.GetConformer()
        lig_conf = self.ligand_no_h.GetConformer()
        lig_positions = lig_conf.GetPositions()
        lig_hyd_indices = sorted(lig_feats["hydrophobes"])

        rows: list[dict] = []
        for res in self.pocket_residues:
            res_key = _residue_key(res)
            if res["residue_name"] not in HYDROPHOBIC_RESIDUES:
                continue
            if res_key in specific_contact_residues:
                continue

            hit_lig_indices: set[int] = set()
            min_dist = float("inf")
            min_prot_name = ""
            min_prot_pos = None

            for p_idx in res["heavy_atom_indices"]:
                p_atom = self.protein.GetAtomWithIdx(p_idx)
                atom_name = p_atom.GetPDBResidueInfo().GetName().strip().upper()
                if atom_name in BACKBONE_ATOMS or p_atom.GetSymbol() not in ("C", "S"):
                    continue
                p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)

                for l_idx in lig_hyd_indices:
                    dist = float(np.linalg.norm(p_pos - lig_positions[l_idx]))
                    if 2.0 <= dist <= hydrophobic_dist:
                        hit_lig_indices.add(l_idx)
                        if dist < min_dist:
                            min_dist = dist
                            min_prot_name = atom_name
                            min_prot_pos = p_pos

            if hit_lig_indices and min_prot_pos is not None:
                sorted_lig = tuple(sorted(hit_lig_indices))
                lig_centroid = lig_positions[list(sorted_lig)].mean(axis=0)
                rows.append(
                    _make_row(
                        res, res_key, False, "Hydrophobic", sorted_lig,
                        min_prot_name, min_dist, None, lig_centroid, min_prot_pos,
                    )
                )

        return rows

    def _detect_clashes(
        self,
        clash_dist: float,
        hbond_pairs: set[tuple[tuple, int]],
    ) -> list[dict]:
        """
        Detect non-H-bonding heavy-atom steric clashes below clash_dist Å.
        """
        prot_conf = self.protein.GetConformer()
        lig_positions = self.ligand_no_h.GetConformer().GetPositions()
        best: dict[tuple[tuple, int], dict] = {}

        for res in self.pocket_residues:
            if res["is_water"]:
                continue
            res_key = _residue_key(res)
            for p_idx in res["heavy_atom_indices"]:
                p_atom = self.protein.GetAtomWithIdx(p_idx)
                atom_name = p_atom.GetPDBResidueInfo().GetName().strip().upper()
                is_bb = atom_name in BACKBONE_ATOMS
                p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)

                for l_idx in range(len(lig_positions)):
                    if (res_key, l_idx) in hbond_pairs:
                        continue
                    dist = float(np.linalg.norm(p_pos - lig_positions[l_idx]))
                    if dist < clash_dist:
                        key = (res_key, l_idx)
                        if key not in best or dist < best[key]["distance"]:
                            best[key] = _make_row(
                                res, res_key, is_bb, "Clash", (l_idx,),
                                atom_name, dist, None, lig_positions[l_idx], p_pos,
                            )

        return list(best.values())

    def _count_prolig_interactions(
        self,
        mol_full: Chem.Mol,
        mol_no_h: Chem.Mol,
        heavy_to_full: dict[int, int],
        pocket_residues: list[dict],
        active_interactions: list[str],
        count_mode: str,
        cutoffs: dict[str, float],
    ) -> dict[tuple[str, tuple], int]:
        """
        Count protein–ligand interactions per (interaction_type, residue_key) for one ligand.
        """
        active_set = set(active_interactions)
        hits: dict[tuple[str, tuple], set] = {}

        def _record(itype: str, res_key: tuple, lig_unit, prot_unit) -> None:
            hits.setdefault((itype, res_key), set()).add(
                lig_unit if count_mode == "ligand" else (lig_unit, prot_unit)
            )

        if not pocket_residues or not active_set:
            return {}

        lig_feats = self._classify_ligand_features(mol_full, mol_no_h, heavy_to_full)
        prot_conf = self.protein.GetConformer()
        lig_conf_full = mol_full.GetConformer()
        lig_conf_heavy = mol_no_h.GetConformer()
        lig_heavy_pos = lig_conf_heavy.GetPositions()

        need_h_or_xb = bool(active_set & {"HBDonor", "HBAcceptor", "XBDonor", "XBAcceptor"})
        prot_donors = []
        prot_acceptors = []
        prot_halogens = []
        if need_h_or_xb:
            for res in pocket_residues:
                res_key = _residue_key(res)
                for p_idx in res["heavy_atom_indices"]:
                    p_atom = self.protein.GetAtomWithIdx(p_idx)
                    anum = p_atom.GetAtomicNum()
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    if anum in (7, 8, 16):
                        atom_name = p_atom.GetPDBResidueInfo().GetName().strip().upper()
                        can_don, can_acc, _ = _protein_hbond_role(
                            res["residue_name"], res["pdb_residue_name"], atom_name, p_atom
                        )
                        if can_don:
                            prot_donors.append((res_key, p_idx, p_atom, p_pos))
                        if can_acc:
                            prot_acceptors.append((res_key, p_idx, p_atom, p_pos))
                    elif "XBAcceptor" in active_set and anum in (17, 35, 53):
                        c_nbrs = [nb for nb in p_atom.GetNeighbors() if nb.GetAtomicNum() == 6]
                        if c_nbrs:
                            c_pos = np.array(prot_conf.GetAtomPosition(c_nbrs[0].GetIdx()), dtype=float)
                            prot_halogens.append((res_key, p_idx, p_pos, c_pos))

        if "HBDonor" in active_set:
            for lig_h_idx, lig_f_idx, _ in lig_feats["donors"]:
                d_atom = mol_full.GetAtomWithIdx(lig_f_idx)
                d_pos = np.array(lig_conf_full.GetAtomPosition(lig_f_idx), dtype=float)
                for res_key, p_idx, a_atom, a_pos in prot_acceptors:
                    dist = float(np.linalg.norm(d_pos - a_pos))
                    if 2.0 <= dist <= cutoffs["hbond_dist"]:
                        ok, _ = _eval_hbond_geometry(
                            d_atom,
                            d_pos,
                            a_pos,
                            lig_conf_full,
                            cutoffs["hbond_angle"],
                            DEFAULT_HBOND_HEAVY_ANGLE,
                            acceptor_atom=a_atom,
                            acceptor_conf=prot_conf,
                        )
                        if ok:
                            _record("HBDonor", res_key, lig_h_idx, p_idx)

        if "HBAcceptor" in active_set:
            for res_key, p_idx, d_atom, d_pos in prot_donors:
                for lig_h_idx, lig_f_idx in lig_feats["acceptors"]:
                    a_atom = mol_full.GetAtomWithIdx(lig_f_idx)
                    a_pos = np.array(lig_conf_full.GetAtomPosition(lig_f_idx), dtype=float)
                    dist = float(np.linalg.norm(d_pos - a_pos))
                    if 2.0 <= dist <= cutoffs["hbond_dist"]:
                        ok, _ = _eval_hbond_geometry(
                            d_atom,
                            d_pos,
                            a_pos,
                            prot_conf,
                            cutoffs["hbond_angle"],
                            DEFAULT_HBOND_HEAVY_ANGLE,
                            acceptor_atom=a_atom,
                            acceptor_conf=lig_conf_full,
                        )
                        if ok:
                            _record("HBAcceptor", res_key, lig_h_idx, p_idx)

        if "XBDonor" in active_set:
            for x_idx, c_idx in lig_feats["halogen_donors"]:
                x_pos = lig_heavy_pos[x_idx]
                c_pos = lig_heavy_pos[c_idx]
                for res_key, p_idx, a_atom, a_pos in prot_acceptors:
                    dist = float(np.linalg.norm(x_pos - a_pos))
                    if 2.0 <= dist <= cutoffs["xbond_dist"] and _eval_xbond_geometry(
                        c_pos,
                        x_pos,
                        a_pos,
                        a_atom,
                        prot_conf,
                        cutoffs["xbond_donor_angle"],
                        cutoffs["xbond_acceptor_angle"],
                    ):
                        _record("XBDonor", res_key, x_idx, p_idx)

        if "XBAcceptor" in active_set and prot_halogens:
            for res_key, p_idx, x_pos, c_pos in prot_halogens:
                for lig_h_idx, _ in lig_feats["acceptors"]:
                    a_atom = mol_no_h.GetAtomWithIdx(lig_h_idx)
                    a_pos = lig_heavy_pos[lig_h_idx]
                    dist = float(np.linalg.norm(x_pos - a_pos))
                    if 2.0 <= dist <= cutoffs["xbond_dist"] and _eval_xbond_geometry(
                        c_pos,
                        x_pos,
                        a_pos,
                        a_atom,
                        lig_conf_heavy,
                        cutoffs["xbond_donor_angle"],
                        cutoffs["xbond_acceptor_angle"],
                    ):
                        _record("XBAcceptor", res_key, lig_h_idx, p_idx)

        need_pi = bool(active_set & {"PiStacking", "FaceToFace", "EdgeToFace", "CationPi"})
        if need_pi:
            for res in pocket_residues:
                if res["residue_name"] not in PROTEIN_AROMATIC_RESIDUES:
                    continue
                res_key = _residue_key(res)
                for p_ring, p_cent, p_norm in self._extract_protein_aromatic_rings(res):
                    if active_set & {"PiStacking", "FaceToFace", "EdgeToFace"}:
                        for l_ring, l_cent, l_norm in lig_feats["arom_rings"]:
                            dist = float(np.linalg.norm(p_cent - l_cent))
                            if dist > cutoffs["pi_stacking_dist"]:
                                continue
                            theta = float(
                                np.degrees(np.arccos(np.clip(abs(float(np.dot(p_norm, l_norm))), 0.0, 1.0)))
                            )
                            is_f2f = dist <= cutoffs["pi_face_dist"] and theta <= DEFAULT_PI_FACE_ANGLE
                            is_e2f = DEFAULT_PI_EDGE_ANGLE <= theta <= 90.0
                            if is_f2f or is_e2f:
                                if "PiStacking" in active_set:
                                    _record("PiStacking", res_key, l_ring, p_ring)
                                if is_f2f and "FaceToFace" in active_set:
                                    _record("FaceToFace", res_key, l_ring, p_ring)
                                if is_e2f and "EdgeToFace" in active_set:
                                    _record("EdgeToFace", res_key, l_ring, p_ring)

                    if "CationPi" in active_set:
                        for l_pos_idx in lig_feats["pos_atoms"]:
                            l_pos = lig_heavy_pos[l_pos_idx]
                            dist = float(np.linalg.norm(p_cent - l_pos))
                            if (
                                dist <= cutoffs["cation_pi_dist"]
                                and _cation_pi_angle(l_pos, p_cent, p_norm) <= cutoffs["cation_pi_angle"]
                            ):
                                _record("CationPi", res_key, l_pos_idx, p_ring)

        if "PiCation" in active_set and lig_feats["arom_rings"]:
            for res in pocket_residues:
                cation_names = _protein_pi_cation_names(res)
                if not cation_names:
                    continue
                res_key = _residue_key(res)
                for p_idx in res["heavy_atom_indices"]:
                    atom_name = self.protein.GetAtomWithIdx(p_idx).GetPDBResidueInfo().GetName().strip().upper()
                    if atom_name not in cation_names:
                        continue
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    for l_ring, l_cent, l_norm in lig_feats["arom_rings"]:
                        dist = float(np.linalg.norm(p_pos - l_cent))
                        if (
                            dist <= cutoffs["cation_pi_dist"]
                            and _cation_pi_angle(p_pos, l_cent, l_norm) <= cutoffs["cation_pi_angle"]
                        ):
                            _record("PiCation", res_key, l_ring, p_idx)

        if active_set & {"Cationic", "Anionic"}:
            for res in pocket_residues:
                pos_names, neg_names = _protein_pos_neg_atom_names(res)
                if not pos_names and not neg_names:
                    continue
                res_key = _residue_key(res)
                for p_idx in res["heavy_atom_indices"]:
                    atom_name = self.protein.GetAtomWithIdx(p_idx).GetPDBResidueInfo().GetName().strip().upper()
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    if "Anionic" in active_set and atom_name in pos_names:
                        for l_idx in lig_feats["neg_atoms"]:
                            if float(np.linalg.norm(p_pos - lig_heavy_pos[l_idx])) <= cutoffs["ionic_dist"]:
                                _record("Anionic", res_key, l_idx, p_idx)
                    if "Cationic" in active_set and atom_name in neg_names:
                        for l_idx in lig_feats["pos_atoms"]:
                            if float(np.linalg.norm(p_pos - lig_heavy_pos[l_idx])) <= cutoffs["ionic_dist"]:
                                _record("Cationic", res_key, l_idx, p_idx)

        if "MetalAcceptor" in active_set and lig_feats["metal_acceptors"]:
            for res in pocket_residues:
                if not res["is_metal"]:
                    continue
                res_key = _residue_key(res)
                for p_idx in res["heavy_atom_indices"]:
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    for l_idx in lig_feats["metal_acceptors"]:
                        if float(np.linalg.norm(p_pos - lig_heavy_pos[l_idx])) <= cutoffs["metal_dist"]:
                            _record("MetalAcceptor", res_key, l_idx, p_idx)

        if "MetalDonor" in active_set and lig_feats["metal_atoms"]:
            for res in pocket_residues:
                if res["is_metal"] or res["is_water"]:
                    continue
                res_key = _residue_key(res)
                for p_idx in res["heavy_atom_indices"]:
                    p_atom = self.protein.GetAtomWithIdx(p_idx)
                    if p_atom.GetAtomicNum() in (7, 8, 16) and p_atom.GetFormalCharge() <= 0:
                        p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                        for l_idx in lig_feats["metal_atoms"]:
                            if float(np.linalg.norm(p_pos - lig_heavy_pos[l_idx])) <= cutoffs["metal_dist"]:
                                _record("MetalDonor", res_key, l_idx, p_idx)

        if "Hydrophobic" in active_set and lig_feats["hydrophobes"]:
            lig_hyd_indices = sorted(lig_feats["hydrophobes"])
            for res in pocket_residues:
                res_key = _residue_key(res)
                for p_idx in res["heavy_atom_indices"]:
                    p_atom = self.protein.GetAtomWithIdx(p_idx)
                    atom_name = p_atom.GetPDBResidueInfo().GetName().strip().upper()
                    if not _is_protein_hydrophobic_atom(res, atom_name, p_atom):
                        continue
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    for l_idx in lig_hyd_indices:
                        dist = float(np.linalg.norm(p_pos - lig_heavy_pos[l_idx]))
                        if 2.0 <= dist <= cutoffs["hydrophobic_dist"]:
                            _record("Hydrophobic", res_key, l_idx, p_idx)

        if "VdWContact" in active_set and len(lig_heavy_pos) > 0:
            pt = Chem.GetPeriodicTable()
            lig_rvdw = np.array(
                [pt.GetRvdw(a.GetAtomicNum()) for a in mol_no_h.GetAtoms()],
                dtype=float,
            )
            tol = cutoffs["vdw_tolerance"]
            for res in pocket_residues:
                res_key = _residue_key(res)
                for p_idx in res["heavy_atom_indices"]:
                    p_atom = self.protein.GetAtomWithIdx(p_idx)
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    p_r = pt.GetRvdw(p_atom.GetAtomicNum())
                    dists = np.linalg.norm(lig_heavy_pos - p_pos, axis=1)
                    for l_idx in np.where(dists <= lig_rvdw + p_r + tol)[0]:
                        _record("VdWContact", res_key, int(l_idx), p_idx)

        if count_mode == "binary":
            return {k: 1 for k, v in hits.items() if v}
        return {k: len(v) for k, v in hits.items() if v}


def _connected_atom_patches(mol: Chem.Mol, atom_indices: list[int]) -> list[tuple[int, ...]]:
    """
    Partition a subset of ligand atom indices into covalently connected components.
    """
    target = set(atom_indices)
    visited: set[int] = set()
    patches: list[tuple[int, ...]] = []

    for start in sorted(target):
        if start in visited:
            continue
        component: list[int] = []
        stack = [start]
        visited.add(start)
        while stack:
            cur = stack.pop()
            component.append(cur)
            for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx in target and nidx not in visited:
                    visited.add(nidx)
                    stack.append(nidx)
        patches.append(tuple(sorted(component)))

    return patches


def _residue_key(res_entry: dict) -> tuple[str, int, str, str]:
    """Return a unique hashable key for a pocket residue entry."""
    return (
        res_entry["chain"],
        res_entry["residue_number"],
        res_entry["insertion_code"],
        res_entry["pdb_residue_name"],
    )


def _make_row(
    res_entry: dict,
    res_key: tuple,
    is_backbone: bool,
    interaction_type: str,
    ligand_atom_indices: tuple[int, ...],
    protein_atom_name: str,
    distance: float,
    angle: Optional[float],
    ligand_coords: np.ndarray,
    protein_coords: np.ndarray,
) -> dict:
    """Assemble one standardized interaction record."""
    return {
        "_res_key": res_key,
        "residue_name": res_entry["residue_name"],
        "residue_number": res_entry["residue_number"],
        "chain": res_entry["chain"],
        "is_backbone": bool(is_backbone),
        "interaction_type": interaction_type,
        "ligand_atom_indices": tuple(ligand_atom_indices),
        "protein_atom_name": protein_atom_name,
        "distance": round(float(distance), 3),
        "angle": round(float(angle), 2) if angle is not None else None,
        "ligand_coords": tuple(round(float(v), 4) for v in ligand_coords),
        "protein_coords": tuple(round(float(v), 4) for v in protein_coords),
    }


def _ring_normal(coords: np.ndarray) -> np.ndarray:
    """Compute a unit normal vector for a planar ring via SVD."""
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered)
    normal = vh[-1]
    norm = float(np.linalg.norm(normal))
    return normal / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])


def _angle_deg(a: np.ndarray, vertex: np.ndarray, b: np.ndarray) -> float:
    """Compute the angle A-vertex-B in degrees."""
    v1 = a - vertex
    v2 = b - vertex
    denom = float(np.linalg.norm(v1) * np.linalg.norm(v2))
    if denom < 1e-12:
        return 0.0
    cos_ang = float(np.clip(np.dot(v1, v2) / denom, -1.0, 1.0))
    return float(np.degrees(np.arccos(cos_ang)))


def _eval_hbond_geometry(
    donor_atom: Chem.Atom,
    donor_pos: np.ndarray,
    acceptor_pos: np.ndarray,
    donor_conf: Chem.Conformer,
    min_h_angle: float,
    min_heavy_angle: float,
    acceptor_atom: Optional[Chem.Atom] = None,
    acceptor_conf: Optional[Chem.Conformer] = None,
    min_acceptor_angle: float = DEFAULT_HBOND_ACCEPTOR_ANGLE,
    max_h_dist: float = DEFAULT_HBOND_H_DIST,
) -> tuple[bool, Optional[float]]:
    """
    Validate hydrogen-bond angle using explicit H neighbors (and optional acceptor heavy-neighbor angle)
    when present, falling back to heavy-atom antecedent geometry on unprotonated structures.
    """
    h_neighbors = [nbr for nbr in donor_atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
    if h_neighbors:
        acc_heavy_pos = (
            [
                np.array(acceptor_conf.GetAtomPosition(nbr.GetIdx()), dtype=float)
                for nbr in acceptor_atom.GetNeighbors()
                if nbr.GetAtomicNum() > 1
            ]
            if acceptor_atom is not None and acceptor_conf is not None
            else []
        )
        best_ang = -1.0
        best_pass_ang = -1.0
        best_ok = False
        for h_atom in h_neighbors:
            h_pos = np.array(donor_conf.GetAtomPosition(h_atom.GetIdx()), dtype=float)
            h_acc_dist = float(np.linalg.norm(h_pos - acceptor_pos))
            ang = _angle_deg(donor_pos, h_pos, acceptor_pos)
            acc_ok = True
            if acc_heavy_pos:
                acc_angles = [_angle_deg(h_pos, acceptor_pos, p) for p in acc_heavy_pos]
                acc_ok = float(np.mean(acc_angles)) >= min_acceptor_angle and min(acc_angles) >= 75.0
            ok_h = h_acc_dist <= max_h_dist and ang >= min_h_angle and acc_ok
            if ang > best_ang:
                best_ang = ang
            if ok_h and ang > best_pass_ang:
                best_pass_ang = ang
                best_ok = True
        return best_ok, (best_pass_ang if best_ok else best_ang)

    heavy_neighbors = [nbr for nbr in donor_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if not heavy_neighbors:
        return True, None

    angles = [
        _angle_deg(
            np.array(donor_conf.GetAtomPosition(nbr.GetIdx()), dtype=float),
            donor_pos,
            acceptor_pos,
        )
        for nbr in heavy_neighbors
    ]
    mean_ang = float(np.mean(angles))
    passed = mean_ang >= min_heavy_angle and min(angles) >= 75.0
    return passed, mean_ang


def _cation_pi_angle(cation_pos: np.ndarray, ring_cent: np.ndarray, ring_norm: np.ndarray) -> float:
    """Return the angle (0..90 degrees) between the ring normal and the centroid-to-cation vector."""
    vec = cation_pos - ring_cent
    dist = float(np.linalg.norm(vec))
    if dist < 1e-8:
        return 0.0
    cos_ang = abs(float(np.dot(ring_norm, vec / dist)))
    return float(np.degrees(np.arccos(np.clip(cos_ang, 0.0, 1.0))))


def _eval_xbond_geometry(
    c_pos: np.ndarray,
    x_pos: np.ndarray,
    a_pos: np.ndarray,
    acceptor_atom: Optional[Chem.Atom],
    acceptor_conf: Optional[Chem.Conformer],
    min_donor_angle: float,
    min_acceptor_angle: float,
) -> bool:
    """Validate halogen-bond C-X...A and optional X...A-R angles."""
    if _angle_deg(c_pos, x_pos, a_pos) < min_donor_angle:
        return False
    if acceptor_atom is not None and acceptor_conf is not None:
        acc_heavy = [nbr for nbr in acceptor_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if acc_heavy:
            acc_angles = [
                _angle_deg(x_pos, a_pos, np.array(acceptor_conf.GetAtomPosition(nbr.GetIdx()), dtype=float))
                for nbr in acc_heavy
            ]
            if float(np.mean(acc_angles)) < min_acceptor_angle or min(acc_angles) < 75.0:
                return False
    return True


def _protein_pos_neg_atom_names(res: dict) -> tuple[set[str], set[str]]:
    """Return sets of positively and negatively charged PDB atom names for a residue."""
    canon = res["residue_name"]
    pdb_res = res["pdb_residue_name"]
    pos_names = set(PROTEIN_POS_ATOMS.get(canon, ()))
    if pdb_res == "HIP":
        pos_names.update({"ND1", "NE2"})
    neg_names = set(PROTEIN_NEG_ATOMS.get(canon, ()))
    if pdb_res in ("ASH", "GLH"):
        neg_names.clear()
    return pos_names, neg_names


def _protein_pi_cation_names(res: dict) -> set[str]:
    """Return PDB atom names of protein cationic centers capable of PiCation interactions."""
    canon = res["residue_name"]
    if canon == "LYS":
        return {"NZ"}
    if canon == "ARG":
        return {"CZ"}
    if res["pdb_residue_name"] == "HIP":
        return {"ND1", "NE2"}
    return set()


def _is_protein_hydrophobic_atom(res: dict, atom_name: str, p_atom: Chem.Atom) -> bool:
    """Return True if a protein heavy atom is a nonpolar sidechain carbon or MET/CYS sulfur."""
    if res["is_water"] or res["is_metal"] or atom_name in BACKBONE_ATOMS:
        return False
    anum = p_atom.GetAtomicNum()
    if anum == 6:
        return not any(nb.GetAtomicNum() in (7, 8, 15) for nb in p_atom.GetNeighbors())
    if anum == 16:
        return res["residue_name"] in ("MET", "CYS") and res["pdb_residue_name"] != "CYM"
    return False


_RES_TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*?)(-?\d+)([A-Za-z]*)$")


def _parse_residue_tag_sort_key(res_tag: str) -> tuple:
    """Parse '<RES><NUM>[CHAIN]' into a sequence-ordered sort key."""
    match = _RES_TAG_RE.match(str(res_tag))
    if not match:
        return (float("inf"), "", str(res_tag), str(res_tag))
    res_name, res_num_str, suffix = match.groups()
    return (int(res_num_str), suffix, res_name, str(res_tag))


def _resolve_color_palette(base_palette: dict, color: Optional[dict]) -> dict[str, tuple[float, float, float]]:
    """Copy base_palette and apply optional RGB tuple or hex/named color overrides."""
    palette = {k: tuple(float(c) for c in v[:3]) for k, v in base_palette.items()}
    if color is None:
        return palette
    if not isinstance(color, dict):
        raise ValueError("color must be a dictionary mapping interaction types to colors.")
    for key, val in color.items():
        if isinstance(val, (tuple, list)) and len(val) in (3, 4):
            palette[key] = tuple(float(c) for c in val[:3])
        else:
            palette[key] = tuple(color_convert(val))
    return palette


def _finalize_matplotlib_figure(
    fig,
    savepath: Optional[str],
    dpi: int,
    show: bool,
):
    """Save a Matplotlib figure when savepath is set and either display or return it."""
    if savepath is not None:
        savepath_str = str(savepath).strip()
        if not savepath_str:
            raise ValueError("savepath must be a non-empty file path when provided.")
        parent_dir = os.path.dirname(savepath_str)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        fig.savefig(savepath_str, dpi=dpi, bbox_inches="tight", pad_inches=0.1, facecolor="white")

    if show:
        plt.show()
        return None
    return fig


def _rgb_to_hex(rgb: tuple[float, float, float]) -> str:
    """Convert normalized (0..1) RGB tuple into a lowercase CSS hex color string."""
    return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(255 * c)))) for c in rgb[:3])


def _prepare_view_payload(
    pli: "PLI",
    pocket_repr: str,
    protein_repr: Optional[str],
    show_labels: bool,
    show_distances: bool,
    split_hbonds: bool,
    include_pharmacophore: bool,
    color: Optional[dict],
    pdb_data: Optional[str] = None,
    include_pdb: bool = True,
    name: Optional[str] = None,
) -> dict:
    """Build backend-neutral 3D receptor/ligand/interaction geometry specs for PLI.view()."""
    palette = _resolve_color_palette(PLI_COLORS, color)
    if include_pdb and pdb_data is None:
        with open(pli.protein_pdb, "r") as f:
            pdb_data = f.read()
    lig_molblock = Chem.MolToMolBlock(pli.ligand_no_h)

    prot_conf = pli.protein.GetConformer()
    res_lookup = {
        (res["chain"], int(res["residue_number"])): res
        for res in pli.pocket_residues
    }

    chains = {rec["chain"] for rec in pli._interaction_records if rec.get("chain")}
    multi_chain = len(chains) > 1

    seen_res: dict[tuple[str, int], dict] = {}
    for rec in pli._interaction_records:
        rk = (rec["chain"], int(rec["residue_number"]))
        if rk not in seen_res:
            seen_res[rk] = rec

    pocket_sels: list[dict] = []
    residue_labels: list[dict] = []
    for (ch, rnum) in sorted(seen_res.keys(), key=lambda x: (x[0], x[1])):
        rec = seen_res[(ch, rnum)]
        sel = {"model": 0, "resi": rnum}
        if ch:
            sel["chain"] = ch
        pocket_sels.append(sel)

        if show_labels:
            res_entry = res_lookup.get((ch, rnum))
            anchor = None
            if res_entry is not None:
                for p_idx in res_entry["heavy_atom_indices"]:
                    info = pli.protein.GetAtomWithIdx(p_idx).GetPDBResidueInfo()
                    if info is not None and info.GetName().strip().upper() == "CA":
                        pos = prot_conf.GetAtomPosition(p_idx)
                        anchor = (float(pos.x), float(pos.y), float(pos.z))
                        break
                if anchor is None and res_entry["heavy_atom_indices"]:
                    pts = np.mean([prot_conf.GetAtomPosition(i) for i in res_entry["heavy_atom_indices"]], axis=0)
                    anchor = (float(pts[0]), float(pts[1]), float(pts[2]))
            if anchor is None:
                anchor = tuple(float(c) for c in rec["protein_coords"])

            residue_labels.append({
                "text": _format_poseview_label(
                    rec["residue_name"], rnum, ch, include_chain=multi_chain
                ),
                "x": anchor[0],
                "y": anchor[1],
                "z": anchor[2],
            })

    dash_types = {"HBond_Donor", "HBond_Acceptor", "Pi_Stacking", "Cation_Pi", "Ionic", "Clash"}
    cylinders: list[dict] = []
    ring_spheres: list[dict] = []
    distance_labels: list[dict] = []

    for rec in pli._interaction_records:
        itype = rec["interaction_type"]
        if itype not in dash_types:
            continue
        color_key = (
            (itype if split_hbonds else "HBond")
            if itype in ("HBond_Donor", "HBond_Acceptor")
            else itype
        )
        hex_col = _rgb_to_hex(palette.get(color_key, (0.2, 0.2, 0.2)))
        lx, ly, lz = [float(v) for v in rec["ligand_coords"]]
        px, py, pz = [float(v) for v in rec["protein_coords"]]
        cylinders.append({
            "start": {"x": lx, "y": ly, "z": lz},
            "end": {"x": px, "y": py, "z": pz},
            "radius": 0.06,
            "dashed": True,
            "dashLength": 0.25,
            "gapLength": 0.18,
            "color": hex_col,
        })
        if len(rec["ligand_atom_indices"]) > 1:
            ring_spheres.append({"x": lx, "y": ly, "z": lz, "radius": 0.18, "color": hex_col})
        if rec.get("protein_atom_name") == "RING":
            ring_spheres.append({"x": px, "y": py, "z": pz, "radius": 0.18, "color": hex_col})
        if show_distances and rec.get("distance") is not None:
            distance_labels.append({
                "text": f"{float(rec['distance']):.1f} Å",
                "x": 0.5 * (lx + px),
                "y": 0.5 * (ly + py),
                "z": 0.5 * (lz + pz),
                "color": hex_col,
            })

    pharm_feats: list[dict] = []
    if include_pharmacophore:
        from pharmacophore.draw import _build_view_feature_specs

        feats = pli._to_pharmacophore()
        pharm_feats = _build_view_feature_specs(feats, INTERACTIVE_COLORS)

    payload = {
        "lig_molblock": lig_molblock,
        "protein_repr": protein_repr,
        "pocket_repr": pocket_repr,
        "pocket_sels": pocket_sels,
        "residue_labels": residue_labels,
        "cylinders": cylinders,
        "ring_spheres": ring_spheres,
        "distance_labels": distance_labels,
        "pharm_feats": pharm_feats,
    }
    if include_pdb:
        payload["pdb_data"] = pdb_data
    if name is not None:
        payload["name"] = name
    return payload


def _render_pli_jupyter(payload: dict, window: tuple[int, int]):
    """Render an interactive py3Dmol widget in Jupyter for a PLI complex."""
    viewer = py3Dmol.view(width=window[0], height=window[1])
    viewer.setBackgroundColor("white")
    viewer.addModel(payload["pdb_data"], "pdb")
    if payload["protein_repr"] == "cartoon":
        viewer.setStyle({"model": 0}, {"cartoon": {"color": "#CBD5E1", "opacity": 0.65}})
    elif payload["protein_repr"] == "line":
        viewer.setStyle({"model": 0}, {"line": {"color": "#CBD5E1", "opacity": 0.45}})
    else:
        viewer.setStyle({"model": 0}, {})

    pocket_style = (
        {"stick": {"colorscheme": "grayCarbon", "radius": 0.16}}
        if payload["pocket_repr"] == "stick"
        else {"line": {"colorscheme": "grayCarbon", "linewidth": 2.0}}
    )
    for sel in payload["pocket_sels"]:
        viewer.addStyle(sel, pocket_style)

    viewer.addModel(payload["lig_molblock"], "mol")
    viewer.setStyle({"model": 1}, {"stick": {"colorscheme": "cyanCarbon", "radius": 0.20}})

    for cyl in payload["cylinders"]:
        viewer.addCylinder(cyl)

    for sp in payload["ring_spheres"]:
        viewer.addSphere({
            "center": {"x": sp["x"], "y": sp["y"], "z": sp["z"]},
            "radius": sp["radius"],
            "color": sp["color"],
            "opacity": 0.9,
        })

    for lbl in payload["residue_labels"]:
        viewer.addLabel(lbl["text"], {
            "position": {"x": lbl["x"], "y": lbl["y"], "z": lbl["z"]},
            "fontSize": 11,
            "fontColor": "#1E293B",
            "showBackground": True,
            "backgroundColor": "#F8FAFC",
            "backgroundOpacity": 0.85,
        })

    for dl in payload["distance_labels"]:
        viewer.addLabel(dl["text"], {
            "position": {"x": dl["x"], "y": dl["y"], "z": dl["z"]},
            "fontSize": 10,
            "fontColor": "white",
            "showBackground": True,
            "backgroundColor": dl["color"],
            "backgroundOpacity": 0.9,
        })

    if payload["pharm_feats"]:
        from pharmacophore.draw import _build_sphere_wedge_mesh

        for feat in payload["pharm_feats"]:
            if feat["phi0"] is None:
                viewer.addSphere({
                    "center": {"x": feat["x"], "y": feat["y"], "z": feat["z"]},
                    "radius": feat["radius"],
                    "color": feat["color"],
                    "opacity": 0.35,
                })
            else:
                mesh = _build_sphere_wedge_mesh(
                    feat["x"], feat["y"], feat["z"], feat["radius"], feat["phi0"], feat["phi1"]
                )
                viewer.addCustom({**mesh, "color": feat["color"], "opacity": 0.35})

    viewer.zoomTo({"model": 1})
    return viewer.show()


def _render_pli_marimo(
    pdb_data: str,
    lig_payloads: list[dict],
    protein_repr: Optional[str],
    pocket_repr: str,
    window: tuple[int, int],
    prefix: str,
):
    """Render a self-contained 3Dmol.js Marimo iframe for one or more PLI complexes."""
    try:
        import marimo as mo
    except Exception as exc:
        return exc

    width, height = window
    sel_style = ' style="display:none"' if len(lig_payloads) <= 1 else ""
    html = f"""
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <div style="font:13px sans-serif;display:flex;gap:8px;align-items:center;margin-bottom:8px">
      <select id="sel"{sel_style}></select>
      <select id="scale" title="Export resolution">
        <option value="1">1x</option><option value="2" selected>2x</option><option value="4">4x</option>
      </select>
      <label><input type="checkbox" id="transparent"> transparent</label>
      <span style="display:inline-flex;gap:2px">
        <button data-pan="-20,0">&larr;</button><button data-pan="20,0">&rarr;</button>
        <button data-pan="0,-20">&uarr;</button><button data-pan="0,20">&darr;</button>
        <button id="recenter" title="Recenter">&#8982;</button>
      </span>
      <button id="shot" title="Save PNG" style="display:flex;align-items:center;justify-content:center;
              width:28px;height:28px;padding:0;cursor:pointer;border:1px solid #ccc;border-radius:4px;background:#fff">
        <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor"
             stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
          <path d="M19 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h11l5 5v11a2 2 0 0 1-2 2z"/>
          <polyline points="17 21 17 13 7 13 7 21"/><polyline points="7 3 7 8 15 8"/>
        </svg>
      </button>
    </div>
    <div id="viewer" style="width:{width}px;height:{height}px;position:relative"></div>
    <script>
    const PDB_DATA = {json.dumps(pdb_data)};
    const LIGANDS = {json.dumps(lig_payloads)};
    const PROTEIN_REPR = {json.dumps(protein_repr)};
    const POCKET_REPR = {json.dumps(pocket_repr)};
    const box = document.getElementById("viewer");
    const sel = document.getElementById("sel");
    LIGANDS.forEach((d, i) => sel.add(new Option(d.name || ("Molecule " + (i + 1)), i)));
    const viewer = $3Dmol.createViewer(box, {{backgroundColor: "white"}});

    function addWedge(viewer, f) {{
      const nLat = 14, nLon = 12, stride = nLon + 1;
      const verts = [], norms = [], faces = [];
      for (let i = 0; i <= nLat; i++) {{
        const th = Math.PI * i / nLat, st = Math.sin(th), ct = Math.cos(th);
        for (let j = 0; j <= nLon; j++) {{
          const ph = f.phi0 + (f.phi1 - f.phi0) * j / nLon;
          const nx = st * Math.sin(ph), ny = ct, nz = st * Math.cos(ph);
          verts.push(new $3Dmol.Vector3(f.x + f.radius * nx, f.y + f.radius * ny, f.z + f.radius * nz));
          norms.push(new $3Dmol.Vector3(nx, ny, nz));
        }}
      }}
      for (let i = 0; i < nLat; i++) {{
        for (let j = 0; j < nLon; j++) {{
          const v00 = i * stride + j, v10 = (i + 1) * stride + j;
          const v11 = (i + 1) * stride + (j + 1), v01 = i * stride + (j + 1);
          if (i > 0) faces.push(v00, v11, v01);
          if (i < nLat - 1) faces.push(v00, v10, v11);
        }}
      }}
      viewer.addCustom({{vertexArr: verts, normalArr: norms, faceArr: faces, color: f.color, opacity: 0.35}});
    }}

    function render(i) {{
      const d = LIGANDS[i];
      viewer.removeAllModels();
      viewer.removeAllShapes();
      viewer.removeAllLabels();

      viewer.addModel(PDB_DATA, "pdb");
      if (PROTEIN_REPR === "cartoon") {{
        viewer.setStyle({{model: 0}}, {{cartoon: {{color: "#CBD5E1", opacity: 0.65}}}});
      }} else if (PROTEIN_REPR === "line") {{
        viewer.setStyle({{model: 0}}, {{line: {{color: "#CBD5E1", opacity: 0.45}}}});
      }} else {{
        viewer.setStyle({{model: 0}}, {{}});
      }}

      const pocketStyle = POCKET_REPR === "stick"
        ? {{stick: {{colorscheme: "grayCarbon", radius: 0.16}}}}
        : {{line: {{colorscheme: "grayCarbon", linewidth: 2.0}}}};
      d.pocket_sels.forEach(s => viewer.addStyle(s, pocketStyle));

      viewer.addModel(d.lig_molblock, "mol");
      viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: "cyanCarbon", radius: 0.20}}}});

      d.cylinders.forEach(c => viewer.addCylinder(c));
      d.ring_spheres.forEach(sp => {{
        viewer.addSphere({{center: {{x: sp.x, y: sp.y, z: sp.z}}, radius: sp.radius, color: sp.color, opacity: 0.9}});
      }});
      d.residue_labels.forEach(lbl => {{
        viewer.addLabel(lbl.text, {{
          position: {{x: lbl.x, y: lbl.y, z: lbl.z}},
          fontSize: 11, fontColor: "#1E293B", showBackground: true,
          backgroundColor: "#F8FAFC", backgroundOpacity: 0.85
        }});
      }});
      d.distance_labels.forEach(dl => {{
        viewer.addLabel(dl.text, {{
          position: {{x: dl.x, y: dl.y, z: dl.z}},
          fontSize: 10, fontColor: "white", showBackground: true,
          backgroundColor: dl.color, backgroundOpacity: 0.9
        }});
      }});
      d.pharm_feats.forEach(f => {{
        if (f.phi0 !== null) addWedge(viewer, f);
        else viewer.addSphere({{center: {{x: f.x, y: f.y, z: f.z}}, radius: f.radius, color: f.color, opacity: 0.35}});
      }});

      viewer.zoomTo({{model: 1}});
      viewer.render();
    }}

    sel.onchange = () => render(+sel.value);
    render(0);

    document.querySelectorAll("[data-pan]").forEach(b => b.onclick = () => {{
      const [dx, dy] = b.dataset.pan.split(",").map(Number);
      viewer.translate(dx, dy, 200);
    }});
    document.getElementById("recenter").onclick = () => viewer.zoomTo({{model: 1}});
    document.getElementById("shot").onclick = () => {{
      const scale = +document.getElementById("scale").value;
      const transparent = document.getElementById("transparent").checked;
      const w = box.style.width, h = box.style.height;
      if (transparent) viewer.setBackgroundColor(0xffffff, 0);
      if (scale !== 1) {{
        box.style.width = parseInt(w) * scale + "px";
        box.style.height = parseInt(h) * scale + "px";
        viewer.resize();
      }}
      viewer.render();
      const uri = viewer.pngURI();
      if (transparent) viewer.setBackgroundColor("white");
      if (scale !== 1) {{ box.style.width = w; box.style.height = h; viewer.resize(); }}
      viewer.render();
      const a = document.createElement("a");
      a.href = uri;
      a.download = LIGANDS.length > 1
        ? `{prefix}_mol_${{+sel.value + 1}}_interactions.png`
        : `{prefix}_interactions.png`;
      document.body.appendChild(a); a.click(); a.remove();
    }};
    </script>"""
    return mo.iframe(html, height=f"{height + 65}px")


def _format_table_residue_tag(res_entry: dict, include_chain: bool) -> str:
    """Format an uppercase 3-letter residue header tag such as 'ASP173' or 'ASP173A'."""
    ch = res_entry["chain"] if (include_chain and res_entry["chain"]) else ""
    ins = res_entry.get("insertion_code", "")
    return f"{res_entry['residue_name']}{res_entry['residue_number']}{ins}{ch}"


def _prepare_barcode_rows(
    df: pd.DataFrame,
    group_by: str,
) -> tuple[list[tuple[str, str, str]], list[float], list[str], list[float]]:
    """Parse, sort, and compute Y-axis tick centers and divider lines for PLI.plot_barcode()."""
    raw_cols = [c for c in df.columns if c not in ("mol_name", "smiles")]
    parsed: list[tuple[str, str, str]] = []
    for col in raw_cols:
        col_str = str(col)
        if "_" in col_str:
            itype, res_tag = col_str.split("_", 1)
        else:
            itype, res_tag = col_str, col_str
        parsed.append((col, itype, res_tag))

    if not parsed:
        return [], [], [], []

    order_map = {name: i for i, name in enumerate(PROLIG_INTERACTION_ORDER)}
    if group_by == "residue":
        parsed.sort(
            key=lambda item: (
                _parse_residue_tag_sort_key(item[2]),
                order_map.get(item[1], len(order_map)),
                item[1],
            )
        )
        ytick_pos: list[float] = []
        ytick_labels: list[str] = []
        dividers: list[float] = []
        start = 0
        while start < len(parsed):
            tag = parsed[start][2]
            end = start + 1
            while end < len(parsed) and parsed[end][2] == tag:
                end += 1
            ytick_pos.append(0.5 * (start + end))
            ytick_labels.append(tag)
            if end < len(parsed):
                dividers.append(float(end))
            start = end
        return parsed, ytick_pos, ytick_labels, dividers

    parsed.sort(
        key=lambda item: (
            order_map.get(item[1], len(order_map)),
            item[1],
            _parse_residue_tag_sort_key(item[2]),
        )
    )
    ytick_pos = [r + 0.5 for r in range(len(parsed))]
    ytick_labels = [item[2] for item in parsed]
    dividers = [float(r) for r in range(1, len(parsed)) if parsed[r][1] != parsed[r - 1][1]]
    return parsed, ytick_pos, ytick_labels, dividers


def _protein_hbond_role(
    canonical_res: str,
    pdb_res: str,
    atom_name: str,
    atom: Chem.Atom,
) -> tuple[bool, bool, bool]:
    """
    Return (can_donate, can_accept, is_backbone) for a protein N/O/S atom.
    """
    if canonical_res in WATER_RESIDUES:
        return True, True, False
    if atom_name == "N":
        return canonical_res != "PRO", False, True
    if atom_name in ("O", "OXT"):
        return False, True, True

    pure_donors = {
        ("LYS", "NZ"),
        ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
        ("TRP", "NE1"),
        ("ASN", "ND2"),
        ("GLN", "NE2"),
    }
    if (canonical_res, atom_name) in pure_donors:
        return True, False, False

    if canonical_res == "ASP" and atom_name in ("OD1", "OD2"):
        return pdb_res == "ASH", True, False
    if canonical_res == "GLU" and atom_name in ("OE1", "OE2"):
        return pdb_res == "GLH", True, False
    if (canonical_res, atom_name) in {("ASN", "OD1"), ("GLN", "OE1"), ("MET", "SD")}:
        return False, True, False

    if (canonical_res, atom_name) in {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH")}:
        return True, True, False
    if canonical_res == "CYS" and atom_name == "SG":
        if pdb_res == "CYX":
            return False, False, False
        if pdb_res == "CYM":
            return False, True, False
        return True, True, False

    if canonical_res == "HIS" and atom_name in ("ND1", "NE2"):
        if any(nbr.GetAtomicNum() == 1 for nbr in atom.GetNeighbors()):
            return True, False, False
        if pdb_res == "HIP":
            return True, False, False
        if pdb_res == "HID":
            return atom_name == "ND1", atom_name == "NE2", False
        if pdb_res == "HIE":
            return atom_name == "NE2", atom_name == "ND1", False
        return True, True, False

    return False, False, False


def _load_protein(protein_pdb: str, add_hs: bool = False) -> Chem.Mol:
    """
    Load a receptor .pdb structure with RDKit while preserving existing hydrogens and handling minor valence quirks.
    """
    if not os.path.isfile(protein_pdb):
        raise FileNotFoundError(f"Protein PDB file not found: {protein_pdb}")
    if os.path.getsize(protein_pdb) == 0:
        raise ValueError(f"Protein PDB file is empty: {protein_pdb}")

    with rdBase.BlockLogs():
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


def _repair_pymol_aromatic_mol(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Repair unsanitized PyMOL-exported aromatic 5-membered ring NH bonds so RDKit can sanitize."""
    repaired = Chem.Mol(mol)
    try:
        repaired.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(repaired)
        for ring in repaired.GetRingInfo().AtomRings():
            if len(ring) != 5:
                continue
            for idx in ring:
                atom = repaired.GetAtomWithIdx(idx)
                if (
                    atom.GetAtomicNum() == 7
                    and atom.GetDegree() == 2
                    and atom.GetTotalNumHs(includeNeighbors=True) == 0
                    and all(b.GetBondType() == Chem.BondType.AROMATIC for b in atom.GetBonds())
                ):
                    atom.SetNumExplicitHs(1)
        Chem.SanitizeMol(repaired)
        return repaired
    except Exception:
        return None


def _prepare_ligand_mol(
    mol: Chem.Mol,
    add_hs: bool = False,
    needs_sanitize: bool = False,
) -> Chem.Mol:
    """Validate 3D conformer presence, repair unsanitized PyMOL aromatic SDFs if needed, and optionally add explicit hydrogens."""
    if mol.GetNumConformers() == 0:
        raise ValueError(
            "Ligand has no 3D conformer. Supply a 3D bound/docked SDF or embed 3D coordinates first."
        )

    prop_name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    if needs_sanitize:
        repaired = _repair_pymol_aromatic_mol(mol)
        if repaired is None:
            raise ValueError(
                "RDKit could not sanitize the SDF ligand. Ensure the SDF has valid bond orders and valence."
            )
        mol = repaired
    if add_hs:
        mol = Chem.AddHs(mol, addCoords=True)

    if prop_name and not (mol.HasProp("_Name") and mol.GetProp("_Name")):
        mol.SetProp("_Name", prop_name)
    return mol


def _read_sdf_records(sdf_path: str, add_hs: bool = False) -> list[Chem.Mol]:
    """Read and prepare all SDF records, falling back to sanitize=False repair when needed."""
    with rdBase.BlockLogs():
        supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
        num_records = len(supplier)

    if num_records == 0:
        raise ValueError(f"No molecule records found in SDF file: {sdf_path}")

    raw_supplier = None
    mols: list[Chem.Mol] = []
    for idx in range(num_records):
        needs_sanitize = False
        with rdBase.BlockLogs():
            mol = supplier[idx]
        if mol is None:
            if raw_supplier is None:
                raw_supplier = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
            mol = raw_supplier[idx]
            needs_sanitize = True
        if mol is None:
            raise ValueError(f"RDKit failed to parse ligand at index {idx} in {sdf_path}.")
        mols.append(
            _prepare_ligand_mol(
                mol,
                add_hs=add_hs,
                needs_sanitize=needs_sanitize,
            )
        )
    return mols


def _load_ligands(
    ligand: Union[str, os.PathLike, Chem.Mol, list[Chem.Mol], tuple[Chem.Mol, ...]],
    add_hs: bool = False,
) -> list[Chem.Mol]:
    """Load one or more 3D ligands from an .sdf file, an RDKit Chem.Mol, or a sequence of Chem.Mol."""
    if isinstance(ligand, Chem.Mol):
        return [_prepare_ligand_mol(Chem.Mol(ligand), add_hs=add_hs)]
    if isinstance(ligand, (str, os.PathLike)):
        sdf_path = str(ligand)
        if not os.path.isfile(sdf_path):
            raise FileNotFoundError(f"Ligand SDF file not found: {sdf_path}")
        return _read_sdf_records(sdf_path, add_hs=add_hs)

    try:
        items = list(ligand)
    except TypeError as exc:
        raise TypeError(
            "ligand must be a path to an .sdf file, an RDKit Chem.Mol, or a sequence of Chem.Mol."
        ) from exc

    if not items:
        raise ValueError("ligand sequence cannot be empty.")
    out: list[Chem.Mol] = []
    for idx, item in enumerate(items):
        if not isinstance(item, Chem.Mol):
            raise TypeError(f"Expected RDKit Chem.Mol at ligand[{idx}], got {type(item).__name__}.")
        out.append(_prepare_ligand_mol(Chem.Mol(item), add_hs=add_hs))
    return out


def _dropdown_labels(names: list[str]) -> list[str]:
    """Build dropdown option labels, prefixing duplicate molecule names with 0-based indices."""
    counts: dict[str, int] = {}
    for name in names:
        counts[name] = counts.get(name, 0) + 1
    return [name if counts[name] == 1 else f"{i}: {name}" for i, name in enumerate(names)]


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


def _format_poseview_label(
    residue_name: str,
    residue_number: int,
    chain: str,
    include_chain: bool = True,
) -> str:
    """
    Format a 2D protein-ligand interaction figure with residue labels such as 'Phe97A' or 'HOH12A'.
    """
    ch = chain if (include_chain and chain) else ""
    canon = AMINO_ACID_ALIASES.get(residue_name.upper(), residue_name.upper())
    if canon in STANDARD_AMINO_ACIDS:
        return f"{canon.capitalize()}{residue_number}{ch}"
    return f"{residue_name.upper()}{residue_number}{ch}"


def _build_residue_2d_template(
    residue_name: str,
    res_label: str,
    is_backbone: bool,
    interaction_type: str,
    protein_atom_name: str,
    target_bond_len: float,
) -> dict:
    """
    Build a 2D RDKit mini-structure template representing an interacting amino acid residue, water, or metal.
    """
    canon = AMINO_ACID_ALIASES.get(residue_name.upper(), residue_name.upper())
    is_prot_donor = interaction_type == "HBond_Acceptor"
    is_arom = interaction_type == "Pi_Stacking" or (
        interaction_type == "Cation_Pi" and protein_atom_name == "RING"
    )

    contact_idx: Optional[int] = None
    anchor_heavy_idx: Optional[int] = None
    stub_idx = 0
    ring_indices: Optional[tuple[int, ...]] = None
    label_map: dict[int, str] = {}
    wedge_bond: Optional[tuple[int, int]] = None
    aux_contacts: dict[str, int] = {}

    if canon in WATER_RESIDUES:
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx, stub_idx = "[H]O[H]", 0, 1, 2
        else:
            smi, contact_idx, anchor_heavy_idx, stub_idx = "O", 0, 0, 0
    elif canon in METAL_RESIDUES:
        smi, contact_idx, anchor_heavy_idx, stub_idx = "[*]", 0, 0, 0
        label_map = {0: res_label}
    elif is_backbone:
        if is_prot_donor or protein_atom_name == "N":
            smi = "[*]C(=O)[CH]([*])N([*])[H]"
            contact_idx, anchor_heavy_idx, stub_idx = 7, 5, 4
            label_map = {0: "R", 4: res_label, 6: "R"}
            wedge_bond = (3, 4)
            aux_contacts = {"N": 7, "O": 2}
        else:
            smi = "[*]C(=O)[CH]([*])N[*]"
            contact_idx, anchor_heavy_idx, stub_idx = 2, 1, 4
            label_map = {0: "R", 4: res_label, 6: "R"}
            wedge_bond = (3, 4)
            aux_contacts = {"O": 2, "N": 5}
    elif is_arom:
        arom_specs = {
            "PHE": ("[*]c1ccccc1", (1, 2, 3, 4, 5, 6)),
            "TYR": ("[*]c1ccc(O)cc1", (1, 2, 3, 4, 6, 7)),
            "HIS": ("[*]c1c[nH]cn1", (1, 2, 3, 4, 5)),
            "TRP": ("[*]c1c[nH]c2ccccc12", (4, 5, 6, 7, 8, 9)),
        }
        smi, ring_indices = arom_specs.get(canon, arom_specs["PHE"])
        stub_idx = 0
        label_map = {0: res_label}
    elif canon == "LYS":
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx = "[*]C[NH2+][H]", 3, 2
        else:
            smi, contact_idx, anchor_heavy_idx = "[*]C[NH3+]", 2, 1
        label_map = {0: res_label}
    elif canon == "ARG":
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx = "[*]NC(=[NH2+])N[H]", 5, 4
            aux_contacts = {"NE": 1, "NH1": 3, "NH2": 5}
        else:
            smi, contact_idx, anchor_heavy_idx = "[*]NC(=[NH2+])N", 3, 2
            aux_contacts = {"NE": 1, "NH1": 3, "NH2": 4}
        label_map = {0: res_label}
    elif canon in ("ASP", "GLU"):
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx = "[*]C(=O)O[H]", 4, 3
        else:
            smi = "[*]C(=O)[O-]"
            contact_idx = 3 if interaction_type == "Ionic" else 2
            anchor_heavy_idx = 1
            aux_contacts = {"OD1": 2, "OD2": 3, "OE1": 2, "OE2": 3}
        label_map = {0: res_label}
    elif canon in ("ASN", "GLN"):
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx = "[*]C(=O)N([H])[H]", 4, 3
        else:
            smi, contact_idx, anchor_heavy_idx = "[*]C(=O)N", 2, 1
        label_map = {0: res_label}
    elif canon in ("SER", "THR"):
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx = "[*]CO[H]", 3, 2
        else:
            smi, contact_idx, anchor_heavy_idx = "[*]CO", 2, 1
        label_map = {0: res_label}
    elif canon == "CYS":
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx = "[*]CS[H]", 3, 2
        else:
            smi, contact_idx, anchor_heavy_idx = "[*]CS", 2, 1
        label_map = {0: res_label}
    elif canon == "TYR":
        if is_prot_donor:
            smi, contact_idx, anchor_heavy_idx = "[*]c1ccc(O[H])cc1", 6, 5
        else:
            smi, contact_idx, anchor_heavy_idx = "[*]c1ccc(O)cc1", 5, 4
        label_map = {0: res_label}
    elif canon == "HIS":
        if is_prot_donor or interaction_type == "Ionic":
            smi, contact_idx, anchor_heavy_idx = "[*]c1c[n]([H])cn1", 4, 3
        else:
            smi, contact_idx, anchor_heavy_idx = "[*]c1c[nH]cn1", 5, 1
        label_map = {0: res_label}
    elif canon == "TRP":
        smi, contact_idx, anchor_heavy_idx = "[*]c1c[n]([H])c2ccccc12", 4, 3
        label_map = {0: res_label}
    else:
        smi, contact_idx, anchor_heavy_idx = "[*]C", 1, 0
        label_map = {0: res_label}

    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    mol = Chem.MolFromSmiles(smi, ps)
    rdDepictor.Compute2DCoords(mol)

    for idx, text in label_map.items():
        mol.GetAtomWithIdx(idx).SetProp("atomLabel", text)
    if wedge_bond is not None:
        bond = mol.GetBondBetweenAtoms(wedge_bond[0], wedge_bond[1])
        if bond is not None:
            bond.SetBondDir(Chem.BondDir.BEGINWEDGE)

    conf = mol.GetConformer()
    pts = np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())],
        dtype=float,
    )
    blens = [
        float(np.linalg.norm(pts[b.GetBeginAtomIdx()] - pts[b.GetEndAtomIdx()]))
        for b in mol.GetBonds()
    ]
    if blens:
        med = float(np.median(blens))
        if med > 1e-6:
            pts *= target_bond_len / med
            for i in range(mol.GetNumAtoms()):
                conf.SetAtomPosition(i, Point3D(float(pts[i, 0]), float(pts[i, 1]), 0.0))

    label_box_radii = []
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 0:
            lbl = atom.GetProp("atomLabel") if atom.HasProp("atomLabel") else "*"
            if len(lbl) > 1:
                label_box_radii.append((1.45 * target_bond_len, 0.72 * target_bond_len))
            else:
                label_box_radii.append((0.85 * target_bond_len, 0.75 * target_bond_len))
        else:
            label_box_radii.append((0.68 * target_bond_len, 0.68 * target_bond_len))

    return {
        "mol": mol,
        "base_coords": pts,
        "contact_idx": contact_idx,
        "anchor_heavy_idx": anchor_heavy_idx,
        "stub_idx": stub_idx,
        "ring_indices": ring_indices,
        "aux_contacts": aux_contacts,
        "label_box_radii": label_box_radii,
    }


def _resolve_template_contact_idx(tpl: dict, rec: dict) -> int:
    """Resolve the local fragment atom index corresponding to a contact record."""
    if tpl["contact_idx"] is None:
        return 0
    if (
        rec.get("interaction_type") == "Ionic"
        and tpl.get("anchor_heavy_idx") is not None
        and tpl["mol"].GetAtomWithIdx(int(tpl["contact_idx"])).GetAtomicNum() == 1
    ):
        return int(tpl["anchor_heavy_idx"])
    aux = tpl.get("aux_contacts", {})
    pname = rec.get("protein_atom_name", "")
    if pname in aux:
        return aux[pname]
    return int(tpl["contact_idx"])


def _rotate_translate_2d(
    base_pts: np.ndarray,
    origin: np.ndarray,
    src_vec: np.ndarray,
    dst_dir: np.ndarray,
    flip_y: bool,
) -> np.ndarray:
    """Rotate and optionally reflect 2D coordinates around origin so src_vec aligns with dst_dir."""
    sn = float(np.linalg.norm(src_vec))
    dn = float(np.linalg.norm(dst_dir))
    if sn < 1e-8 or dn < 1e-8:
        return base_pts.copy()

    ex = dst_dir / dn
    ey = np.array([-ex[1], ex[0]], dtype=float)

    ax_u = src_vec / sn
    ax_v = np.array([-ax_u[1], ax_u[0]], dtype=float)
    if flip_y:
        ax_v = -ax_v

    shifted = base_pts - origin
    proj_x = shifted @ ax_u
    proj_y = shifted @ ax_v
    return proj_x[:, None] * ex + proj_y[:, None] * ey


def _orient_fragment_coords(
    tpl: dict,
    anchor_xy: np.ndarray,
    theta: float,
    dist: float,
    lig_pts: np.ndarray,
    l_bond: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Place a 2D residue template at (theta, dist) around anchor_xy, pointing its reactive contact atom
    toward the ligand while flipping laterally for optimal ligand clearance.
    """
    u = np.array([np.cos(theta), np.sin(theta)], dtype=float)
    w = -u
    q = anchor_xy + dist * u
    base = tpl["base_coords"]
    mol = tpl["mol"]

    if len(base) == 1:
        return np.tile(q, (1, 1)), q

    best_pts = None
    best_clear = -1e18

    if tpl["ring_indices"] is not None:
        r_idxs = list(tpl["ring_indices"])
        c0 = base[r_idxs].mean(axis=0)
        v_cur = base[tpl["stub_idx"]] - c0
        for flip in (False, True):
            rot = _rotate_translate_2d(base, c0, v_cur, u, flip_y=flip)
            cand = rot + q
            clr = min(
                float(np.linalg.norm(lig_pts - cand[i], axis=1).min())
                for i in range(len(cand))
            )
            if clr > best_clear:
                best_clear = clr
                best_pts = cand
        return best_pts, q

    c_idx = tpl["contact_idx"]
    h_idx = tpl["anchor_heavy_idx"]
    origin = base[h_idx]
    v_cur = base[c_idx] - origin
    is_collinear_h = mol.GetAtomWithIdx(c_idx).GetAtomicNum() == 1 and c_idx != h_idx

    for flip in (False, True):
        rot = _rotate_translate_2d(base, origin, v_cur, w, flip_y=flip)
        if is_collinear_h:
            rot[c_idx] = rot[h_idx] + l_bond * w
        cand = rot + (q - rot[c_idx])
        clr = min(
            float(np.linalg.norm(lig_pts - cand[i], axis=1).min())
            for i in range(len(cand))
        )
        if clr > best_clear:
            best_clear = clr
            best_pts = cand

    return best_pts, q


def _segments_intersect_2d(a0: np.ndarray, a1: np.ndarray, b0: np.ndarray, b1: np.ndarray) -> bool:
    """Return True if two closed 2D line segments a0-a1 and b0-b1 intersect cleanly in their interiors."""
    def _cross(o: np.ndarray, p: np.ndarray, q: np.ndarray) -> float:
        return float((p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0]))

    if (
        max(a0[0], a1[0]) < min(b0[0], b1[0])
        or max(b0[0], b1[0]) < min(a0[0], a1[0])
        or max(a0[1], a1[1]) < min(b0[1], b1[1])
        or max(b0[1], b1[1]) < min(a0[1], a1[1])
    ):
        return False

    c1 = _cross(a0, a1, b0)
    c2 = _cross(a0, a1, b1)
    c3 = _cross(b0, b1, a0)
    c4 = _cross(b0, b1, a1)
    return (c1 * c2 < -1e-6) and (c3 * c4 < -1e-6)


def _point_to_segments_min_dist(points: np.ndarray, segs: list[tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    """Compute the minimum Euclidean distance from each (M, 2) point to a list of 2D line segments."""
    if not segs:
        return np.full(len(points), np.inf, dtype=float)
    min_d2 = np.full(len(points), np.inf, dtype=float)
    for a, b in segs:
        ab = b - a
        ab2 = float(np.dot(ab, ab))
        if ab2 < 1e-12:
            d2 = np.sum((points - a) ** 2, axis=1)
        else:
            t = np.clip(((points - a) @ ab) / ab2, 0.0, 1.0)
            proj = a + t[:, None] * ab
            d2 = np.sum((points - proj) ** 2, axis=1)
        min_d2 = np.minimum(min_d2, d2)
    return np.sqrt(min_d2)


def _relax_2d_fragments(
    frag_specs: list[dict],
    lig_pts: np.ndarray,
    mol_2d: Chem.Mol,
    l_bond: float,
    max_sweeps: int = 80,
) -> None:
    """
    Iteratively relax residue fragment angles and radial distances in 2D molecule units to eliminate
    ligand overlaps, fragment/label collisions, and connector crossings.
    """
    lig_bonds = [
        (lig_pts[b.GetBeginAtomIdx()], lig_pts[b.GetEndAtomIdx()])
        for b in mol_2d.GetBonds()
    ]
    lig_coll_pts = list(lig_pts) + [0.5 * (p0 + p1) for p0, p1 in lig_bonds]
    lig_coll = np.array(lig_coll_pts, dtype=float)

    n_frags = len(frag_specs)
    cache_pts: list[np.ndarray] = []
    cache_qs: list[np.ndarray] = []
    for spec in frag_specs:
        pts, q = _orient_fragment_coords(
            spec["tpl"], spec["anchor_xy"], spec["theta"], spec["d"], lig_pts, l_bond
        )
        cache_pts.append(pts)
        cache_qs.append(q)

    def _layout_energy() -> float:
        energy = 0.0
        # Soft anchor / radial priors
        for spec in frag_specs:
            dtheta = float(np.arctan2(np.sin(spec["theta"] - spec["theta_pref"]), np.cos(spec["theta"] - spec["theta_pref"])))
            dd = (spec["d"] - spec["d_0"]) / l_bond
            energy += 0.25 * (dtheta ** 2) + 0.15 * (dd ** 2)

        # Fragment-vs-ligand collisions & connector-vs-ligand-bond crossings
        for i in range(n_frags):
            spec_i = frag_specs[i]
            pts_i = cache_pts[i]
            radii_i = spec_i["tpl"]["label_box_radii"]
            for a_idx, (ax, ay) in enumerate(pts_i):
                rx, ry = radii_i[a_idx]
                dx = np.abs(lig_coll[:, 0] - ax) / (rx + 0.55 * l_bond)
                dy = np.abs(lig_coll[:, 1] - ay) / (ry + 0.55 * l_bond)
                norm2 = dx ** 2 + dy ** 2
                pen = np.maximum(0.0, 1.0 - np.sqrt(norm2))
                energy += 30.0 * float(np.sum(pen ** 2))

            p_i = spec_i["anchor_xy"]
            q_i = cache_qs[i]
            u_i = q_i - p_i
            p_inset = p_i + 0.2 * u_i
            for b0, b1 in lig_bonds:
                if _segments_intersect_2d(p_inset, q_i, b0, b1):
                    energy += 80.0

        # Fragment-vs-fragment collisions & connector crossings
        for i in range(n_frags):
            pts_i = cache_pts[i]
            rad_i = frag_specs[i]["tpl"]["label_box_radii"]
            p_i, q_i = frag_specs[i]["anchor_xy"], cache_qs[i]
            for j in range(i + 1, n_frags):
                pts_j = cache_pts[j]
                rad_j = frag_specs[j]["tpl"]["label_box_radii"]
                p_j, q_j = frag_specs[j]["anchor_xy"], cache_qs[j]
                if _segments_intersect_2d(p_i, q_i, p_j, q_j):
                    energy += 120.0
                for ai, (xi, yi) in enumerate(pts_i):
                    rxi, ryi = rad_i[ai]
                    for aj, (xj, yj) in enumerate(pts_j):
                        rxj, ryj = rad_j[aj]
                        norm_d = np.hypot((xi - xj) / (rxi + rxj), (yi - yj) / (ryi + ryj))
                        if norm_d < 1.0:
                            energy += 35.0 * float((1.0 - norm_d) ** 2)
                # Connector j passing through fragment i
                if len(pts_i) > 1:
                    d_seg = _point_to_segments_min_dist(pts_i, [(p_j, q_j)])
                    pen_seg = np.maximum(0.0, 0.85 * l_bond - d_seg)
                    energy += 15.0 * float(np.sum(pen_seg ** 2))
                if len(pts_j) > 1:
                    d_seg = _point_to_segments_min_dist(pts_j, [(p_i, q_i)])
                    pen_seg = np.maximum(0.0, 0.85 * l_bond - d_seg)
                    energy += 15.0 * float(np.sum(pen_seg ** 2))

        return energy

    # First check direct connector crossings and flip angles apart if any
    for i in range(n_frags):
        for j in range(i + 1, n_frags):
            if _segments_intersect_2d(
                frag_specs[i]["anchor_xy"], cache_qs[i],
                frag_specs[j]["anchor_xy"], cache_qs[j],
            ):
                ti, tj = frag_specs[i]["theta"], frag_specs[j]["theta"]
                frag_specs[i]["theta"], frag_specs[j]["theta"] = tj, ti
                cache_pts[i], cache_qs[i] = _orient_fragment_coords(
                    frag_specs[i]["tpl"], frag_specs[i]["anchor_xy"], frag_specs[i]["theta"],
                    frag_specs[i]["d"], lig_pts, l_bond,
                )
                cache_pts[j], cache_qs[j] = _orient_fragment_coords(
                    frag_specs[j]["tpl"], frag_specs[j]["anchor_xy"], frag_specs[j]["theta"],
                    frag_specs[j]["d"], lig_pts, l_bond,
                )

    best_energy = _layout_energy()
    if best_energy < 1e-4:
        return

    step_thetas = np.radians([-14.0, -6.0, 0.0, 6.0, 14.0])
    step_ds = np.array([-0.25, 0.0, 0.35], dtype=float) * l_bond

    for _ in range(max_sweeps):
        improved = False
        for k in range(n_frags):
            spec = frag_specs[k]
            cur_th = spec["theta"]
            cur_d = spec["d"]
            min_d = max(1.8 * l_bond, spec["d_0"] - 0.35 * l_bond)
            max_d = spec["d_0"] + 2.4 * l_bond

            best_k = (cur_th, cur_d, cache_pts[k], cache_qs[k], best_energy)
            for dt in step_thetas:
                cand_th = cur_th + float(dt)
                dev = float(np.arctan2(
                    np.sin(cand_th - spec["theta_init"]),
                    np.cos(cand_th - spec["theta_init"]),
                ))
                if abs(dev) > np.radians(78.0):
                    continue
                for dd in step_ds:
                    if dt == 0.0 and dd == 0.0:
                        continue
                    cand_d = float(np.clip(cur_d + float(dd), min_d, max_d))
                    spec["theta"] = cand_th
                    spec["d"] = cand_d
                    cache_pts[k], cache_qs[k] = _orient_fragment_coords(
                        spec["tpl"], spec["anchor_xy"], cand_th, cand_d, lig_pts, l_bond
                    )
                    e = _layout_energy()
                    if e < best_k[4] - 1e-6:
                        best_k = (cand_th, cand_d, cache_pts[k], cache_qs[k], e)
                        improved = True

            spec["theta"], spec["d"], cache_pts[k], cache_qs[k], best_energy = best_k

        if not improved or best_energy < 1e-4:
            break


def _extract_kept_runs(
    poly: np.ndarray,
    keep: np.ndarray,
) -> list[np.ndarray]:
    """Split a contourpy polyline into contiguous kept point runs, merging closed-loop wrap-arounds."""
    n = len(poly)
    if n < 2 or not keep.any():
        return []

    is_closed = n >= 3 and bool(np.allclose(poly[0], poly[-1], atol=1e-3))
    pts = poly[:-1] if is_closed else poly
    mask = keep[:-1] if is_closed else keep

    if mask.all():
        return [poly.copy()]

    runs: list[list[int]] = []
    cur: list[int] = []
    for idx, flag in enumerate(mask):
        if flag:
            cur.append(idx)
        elif cur:
            runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)

    if is_closed and len(runs) >= 2 and runs[0][0] == 0 and runs[-1][-1] == len(pts) - 1:
        merged = runs[-1] + runs[0]
        runs = [merged] + runs[1:-1]

    return [pts[np.array(r)] for r in runs if len(r) >= 3]


def _smooth_polyline(pts: np.ndarray, passes: int = 2) -> np.ndarray:
    """Smooth an open 2D polyline lightly via a 3-point binomial filter while keeping endpoints fixed."""
    out = pts.copy()
    if len(out) < 4:
        return out
    for _ in range(passes):
        nxt = out.copy()
        nxt[1:-1] = 0.25 * out[:-2] + 0.5 * out[1:-1] + 0.25 * out[2:]
        out = nxt
    return out


def _build_hydrophobic_splines_and_labels(
    records: list[dict],
    mol_2d: Chem.Mol,
    n_heavy_lig: int,
    atom_pixels: dict[int, tuple[float, float]],
    frag_atom_pixels: list[tuple[float, float]],
    connectors: list[dict],
    polar_lig_heavy_atoms: set[int],
    bond_px: float,
    canvas_size: tuple[int, int],
    include_chain: bool,
) -> list[dict]:
    """
    Group contacting hydrophobic residues into spatial ligand bays, extract smooth offset contour splines
    via contourpy, and place non-overlapping stacked unboxed residue labels outside each curve.
    """
    hydro_rows = [r for r in records if r["interaction_type"] == "Hydrophobic"]
    if not hydro_rows:
        return []

    row_pts = [
        np.array([atom_pixels[i] for i in r["ligand_atom_indices"]], dtype=float)
        for r in hydro_rows
    ]

    # Connected components grouping hydrophobic rows whose contacted atoms lie within 2.5 * bond_px
    n_rows = len(hydro_rows)
    visited = [False] * n_rows
    components: list[list[int]] = []
    merge_dist = 2.5 * bond_px

    for start in range(n_rows):
        if visited[start]:
            continue
        comp = []
        stack = [start]
        visited[start] = True
        while stack:
            cur = stack.pop()
            comp.append(cur)
            pts_c = row_pts[cur]
            for nxt in range(n_rows):
                if visited[nxt]:
                    continue
                dmin = float(np.linalg.norm(pts_c[:, None, :] - row_pts[nxt][None, :, :], axis=2).min())
                if dmin <= merge_dist:
                    visited[nxt] = True
                    stack.append(nxt)
        components.append(comp)

    width, height = canvas_size
    lig_heavy_pts = np.array([atom_pixels[i] for i in range(n_heavy_lig)], dtype=float)
    lig_segs = [
        (
            np.array(atom_pixels[b.GetBeginAtomIdx()], dtype=float),
            np.array(atom_pixels[b.GetEndAtomIdx()], dtype=float),
        )
        for b in mol_2d.GetBonds()
        if b.GetBeginAtomIdx() < n_heavy_lig and b.GetEndAtomIdx() < n_heavy_lig
    ]
    conn_segs = [
        (np.array(c["lig_point"], dtype=float), np.array(c["prot_point"], dtype=float))
        for c in connectors
    ]
    frag_pts_arr = np.array(frag_atom_pixels, dtype=float) if frag_atom_pixels else np.empty((0, 2))
    lig_all_pts = np.array([atom_pixels[i] for i in range(mol_2d.GetNumAtoms())], dtype=float)

    step = 4.0
    xs = np.arange(0.0, width + step, step, dtype=float)
    ys = np.arange(0.0, height + step, step, dtype=float)
    X, Y = np.meshgrid(xs, ys)
    grid_pts = np.column_stack([X.ravel(), Y.ravel()])

    if lig_segs:
        d_grid = _point_to_segments_min_dist(grid_pts, lig_segs).reshape(X.shape)
    else:
        d_grid = np.linalg.norm(grid_pts[:, None, :] - lig_heavy_pts[None, :, :], axis=2).min(axis=1).reshape(X.shape)

    r_arc = 0.72 * bond_px
    cg = contourpy.contour_generator(x=xs, y=ys, z=d_grid)
    contour_polys = cg.lines(r_arc)

    placed_label_boxes: list[tuple[float, float, float, float]] = []
    patches_out: list[dict] = []

    for comp in components:
        labels: list[str] = []
        seen_lbls: set[str] = set()
        patch_atoms: set[int] = set()
        for r_idx in comp:
            rec = hydro_rows[r_idx]
            patch_atoms.update(rec["ligand_atom_indices"])
            lbl = _format_poseview_label(
                rec["residue_name"],
                rec["residue_number"],
                rec["chain"],
                include_chain=include_chain,
            )
            if lbl not in seen_lbls:
                seen_lbls.add(lbl)
                labels.append(lbl)

        h_pts = np.array([atom_pixels[i] for i in sorted(patch_atoms)], dtype=float)
        other_polar = sorted(polar_lig_heavy_atoms - patch_atoms)
        other_polar_pts = np.array([atom_pixels[i] for i in other_polar], dtype=float) if other_polar else np.empty((0, 2))

        kept_arcs: list[np.ndarray] = []
        for poly in contour_polys:
            d_h = np.linalg.norm(poly[:, None, :] - h_pts[None, :, :], axis=2).min(axis=1)
            keep = d_h <= 1.25 * r_arc
            if len(other_polar_pts):
                d_pol = np.linalg.norm(poly[:, None, :] - other_polar_pts[None, :, :], axis=2).min(axis=1)
                keep &= d_h <= 1.15 * d_pol
            if len(frag_pts_arr):
                d_frag = np.linalg.norm(poly[:, None, :] - frag_pts_arr[None, :, :], axis=2).min(axis=1)
                keep &= d_frag >= 0.78 * bond_px
            if conn_segs:
                d_conn = _point_to_segments_min_dist(poly, conn_segs)
                keep &= d_conn >= 0.55 * bond_px

            for run in _extract_kept_runs(poly, keep):
                arc_len = float(np.sum(np.linalg.norm(np.diff(run, axis=0), axis=1)))
                if arc_len >= 0.55 * bond_px:
                    kept_arcs.append(_smooth_polyline(run, passes=4))

        # Keep up to the 2 longest distinct runs in a patch to keep visual style clean
        kept_arcs.sort(
            key=lambda a: float(np.sum(np.linalg.norm(np.diff(a, axis=0), axis=1))),
            reverse=True,
        )
        kept_arcs = kept_arcs[:2]

        # Fallback outward circular arc if tight crowding removed every contour run
        if not kept_arcs:
            lig_center = lig_heavy_pts.mean(axis=0)
            out_dir = h_pts.mean(axis=0) - lig_center
            nrm = float(np.linalg.norm(out_dir))
            if nrm < 1e-6:
                out_dir = np.array([0.0, 1.0])
            else:
                out_dir /= nrm
            ang0 = float(np.arctan2(out_dir[1], out_dir[0]))
            angs = np.linspace(ang0 - 0.55, ang0 + 0.55, 16)
            most_exp = h_pts[int(np.argmax((h_pts - lig_center) @ out_dir))]
            arc_r = 0.78 * bond_px
            kept_arcs = [np.column_stack([most_exp[0] + arc_r * np.cos(angs), most_exp[1] + arc_r * np.sin(angs)])]

        primary_arc = kept_arcs[0]
        mid_idx = len(primary_arc) // 2
        arc_mid = primary_arc[mid_idx]
        nearest_h = h_pts[int(np.argmin(np.linalg.norm(h_pts - arc_mid, axis=1)))]
        n_vec = arc_mid - nearest_h
        nn = float(np.linalg.norm(n_vec))
        n_out = n_vec / nn if nn > 1e-6 else np.array([0.0, 1.0])

        max_chars = max(len(lbl) for lbl in labels)
        est_w = max_chars * 0.17 * bond_px + 12.0
        est_h = len(labels) * 0.34 * bond_px + 8.0

        if n_out[0] > 0.25:
            ha = "left"
        elif n_out[0] < -0.25:
            ha = "right"
        else:
            ha = "center"

        base_offset = 0.50 * bond_px + 0.5 * est_h * abs(n_out[1]) + 0.15 * est_w * abs(n_out[0])
        best_pos = arc_mid + base_offset * n_out

        # Nudge stacked label outward / sideways if it would overlap fragments, connectors, or earlier labels
        t_vec = np.array([-n_out[1], n_out[0]], dtype=float)
        best_score = -1e18
        for radial_mul in (1.0, 1.25, 1.55, 1.9, 2.3):
            for tang_shift in (0.0, -0.8 * bond_px, 0.8 * bond_px, -1.6 * bond_px, 1.6 * bond_px):
                cand = arc_mid + (base_offset * radial_mul) * n_out + tang_shift * t_vec
                if ha == "left":
                    cx0, cx1 = cand[0], cand[0] + est_w
                elif ha == "right":
                    cx0, cx1 = cand[0] - est_w, cand[0]
                else:
                    cx0, cx1 = cand[0] - 0.5 * est_w, cand[0] + 0.5 * est_w
                cy0, cy1 = cand[1] - 0.5 * est_h, cand[1] + 0.5 * est_h

                pen = 0.0
                # Check ligand + fragment atom overlap against label rectangle
                for obs_pts in (lig_all_pts, frag_pts_arr):
                    if not len(obs_pts):
                        continue
                    ox = np.maximum(0.0, np.maximum(cx0 - obs_pts[:, 0], obs_pts[:, 0] - cx1))
                    oy = np.maximum(0.0, np.maximum(cy0 - obs_pts[:, 1], obs_pts[:, 1] - cy1))
                    min_box_d = float(np.hypot(ox, oy).min())
                    if min_box_d < 0.65 * bond_px:
                        pen += (0.65 * bond_px - min_box_d) * 5.0
                # Check prior label boxes
                for px0, px1, py0, py1 in placed_label_boxes:
                    if not (cx1 < px0 - 6.0 or cx0 > px1 + 6.0 or cy1 < py0 - 6.0 or cy0 > py1 + 6.0):
                        pen += 100.0
                # Check connectors
                if conn_segs:
                    box_samples = np.array(
                        [[cand[0], cand[1]], [cx0, cy0], [cx1, cy0], [cx0, cy1], [cx1, cy1]],
                        dtype=float,
                    )
                    d_c = float(_point_to_segments_min_dist(box_samples, conn_segs).min())
                    if d_c < 0.5 * bond_px:
                        pen += (0.5 * bond_px - d_c) * 4.0

                score = -pen - 0.15 * radial_mul - 0.05 * abs(tang_shift / bond_px)
                if score > best_score:
                    best_score = score
                    best_pos = cand

        if ha == "left":
            placed_label_boxes.append((best_pos[0], best_pos[0] + est_w, best_pos[1] - 0.5 * est_h, best_pos[1] + 0.5 * est_h))
        elif ha == "right":
            placed_label_boxes.append((best_pos[0] - est_w, best_pos[0], best_pos[1] - 0.5 * est_h, best_pos[1] + 0.5 * est_h))
        else:
            placed_label_boxes.append((best_pos[0] - 0.5 * est_w, best_pos[0] + 0.5 * est_w, best_pos[1] - 0.5 * est_h, best_pos[1] + 0.5 * est_h))

        patches_out.append(
            {
                "arcs": kept_arcs,
                "labels": labels,
                "label_pos": (float(best_pos[0]), float(best_pos[1])),
                "label_ha": ha,
                "label_box_size": (est_w, est_h),
            }
        )

    return patches_out


def _atom_glyph_radius_px(atom: Chem.Atom, bond_px: float) -> float:
    """Estimate pixel clearance radius around an atom symbol so dashed connectors stop outside glyphs."""
    anum = atom.GetAtomicNum()
    if anum == 1:
        return 0.28 * bond_px
    if anum == 0:
        lbl = atom.GetProp("atomLabel") if atom.HasProp("atomLabel") else "*"
        return (0.48 if len(lbl) > 2 else 0.32) * bond_px
    if atom.GetTotalNumHs() > 0 or atom.GetFormalCharge() != 0:
        return 0.36 * bond_px
    if anum == 6:
        return 0.14 * bond_px
    return 0.28 * bond_px


def _clip_segment_for_glyphs(
    p0: np.ndarray,
    p1: np.ndarray,
    pad0_px: float,
    pad1_px: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Clip segment endpoints inward along its unit direction without reversing short connectors."""
    vec = p1 - p0
    length = float(np.linalg.norm(vec))
    if length < 1e-6:
        return p0, p1
    u = vec / length
    total_pad = pad0_px + pad1_px
    if total_pad >= 0.8 * length:
        scale = 0.8 * length / max(total_pad, 1e-6)
        pad0_px *= scale
        pad1_px *= scale
    return p0 + pad0_px * u, p1 - pad1_px * u


def _draw_dashed_arrow(
    ax,
    start: np.ndarray,
    end: np.ndarray,
    color,
    bond_px: float,
) -> None:
    """Draw a dashed connector from start to end with a solid directional arrowhead at end."""
    vec = end - start
    seg_len = float(np.linalg.norm(vec))
    if seg_len < 1e-6:
        return
    u = vec / seg_len
    perp = np.array([-u[1], u[0]], dtype=float)

    head_len = min(0.24 * bond_px, 0.35 * seg_len)
    head_half_w = 0.40 * head_len
    shaft_end = end - 0.85 * head_len * u

    ax.add_line(
        Line2D(
            [start[0], shaft_end[0]],
            [start[1], shaft_end[1]],
            linestyle=(0, (4.2, 2.8)),
            linewidth=1.9,
            color=color,
            solid_capstyle="butt",
            zorder=4,
        )
    )
    base_center = end - head_len * u
    tri = np.array(
        [
            end,
            base_center + head_half_w * perp,
            base_center - head_half_w * perp,
        ],
        dtype=float,
    )
    ax.add_patch(
        patches.Polygon(
            tri,
            closed=True,
            facecolor=color,
            edgecolor="none",
            zorder=5,
        )
    )


def _draw_pli_legend(ax, present_types: set[str], palette: dict) -> None:
    """Render a compact 2D protein-ligand interaction legend for present interaction types."""
    legend_order = [
        ("HBond", "H-Bond (D → A)"),
        ("HBond_Donor", "H-Bond (Lig Donor)"),
        ("HBond_Acceptor", "H-Bond (Lig Acceptor)"),
        ("Ionic", "Ionic / Salt Bridge"),
        ("Pi_Stacking", "π–π Stacking"),
        ("Cation_Pi", "Cation–π"),
        ("Hydrophobic", "Hydrophobic"),
        ("Clash", "Steric Clash"),
    ]
    handles = []
    for key, text in legend_order:
        if key not in present_types:
            continue
        col = palette.get(key, (0.2, 0.2, 0.2))
        if key == "Hydrophobic":
            handles.append(Line2D([0], [0], color=col, lw=2.0, linestyle="-", label=text))
        elif key in ("HBond", "HBond_Donor", "HBond_Acceptor"):
            handles.append(
                Line2D([0], [0], color=col, lw=1.8, linestyle=(0, (4.0, 2.8)), marker=">", markersize=5.5, label=text)
            )
        elif key in ("Pi_Stacking", "Cation_Pi"):
            handles.append(
                Line2D([0], [0], color=col, lw=1.8, linestyle=(0, (4.0, 2.8)), marker="o", markersize=5.5, label=text)
            )
        else:
            handles.append(Line2D([0], [0], color=col, lw=1.8, linestyle=(0, (4.0, 2.8)), label=text))

    if handles:
        ax.legend(
            handles=handles,
            loc="upper right",
            fontsize=8.5,
            frameon=True,
            facecolor="white",
            edgecolor="#CBD5E1",
            framealpha=0.95,
        )
