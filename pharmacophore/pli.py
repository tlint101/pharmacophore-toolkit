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
    BACKBONE_ATOMS,
    DEFAULT_CATION_PI_DIST,
    DEFAULT_CLASH_DIST,
    DEFAULT_HBOND_ANGLE,
    DEFAULT_HBOND_DIST,
    DEFAULT_HBOND_HEAVY_ANGLE,
    DEFAULT_HYDROPHOBIC_DIST,
    DEFAULT_IONIC_DIST,
    DEFAULT_PI_STACKING_DIST,
    DEFAULT_POCKET_CUTOFF,
    FEATURE_COLORS,
    HYDROPHOBIC_RESIDUES,
    METAL_RESIDUES,
    PML_COLORS,
    PROTEIN_AROMATIC_RESIDUES,
    PROTEIN_AROMATIC_RING_ATOMS,
    PROTEIN_NEG_ATOMS,
    PROTEIN_POS_ATOMS,
    STANDARD_AMINO_ACIDS,
    WATER_RESIDUES,
    color_convert,
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
        self._interaction_records = []
        self.interactions_df = None
        self.features_list = None

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

    def detect(
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

    def _to_pharmacophore(self, **detect_kwargs) -> list:
        """
        Internal helper to convert detected 3D protein–ligand interactions into pharmacophore-toolkit's
        standard feature list.
        :param detect_kwargs: dict
            Optional keyword arguments forwarded to PLI.detect() (forces recomputation when supplied).
        :return: list
            List of [family: str, atom_indices: tuple, x: float, y: float, z: float] entries.
        """
        if detect_kwargs or self.interactions_df is None:
            self.detect(**detect_kwargs)

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
    ) -> str:
        """
        Write a drag-and-drop PyMOL (.pml) visualization script for the detected protein–ligand interactions.
        :param savepath: str
            Output file path ending in .pml.
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

        if self.interactions_df is None:
            self.detect()

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
                pharm_counts: dict[str, int] = {}
                for fam, _, x, y, z in feats:
                    pharm_counts[fam] = pharm_counts.get(fam, 0) + 1
                    obj_name = f"pli_pharm_{fam}_{pharm_counts[fam]}"
                    f.write(f"pseudoatom {obj_name}, pos=[{x}, {y}, {z}]\n")
                    f.write(f"show spheres, {obj_name}\n")
                    f.write(f"hide nonbonded, {obj_name}\n")
                    f.write(f"color pli_pharm_{fam.lower()}_color, {obj_name}\n")
                    f.write(f"set sphere_scale, {sphere_size}, {obj_name}\n")
                    f.write(f"set sphere_transparency, 0.25, {obj_name}\n\n")
                    group_members.append(obj_name)

            if group_members:
                f.write(f"group pli_interactions, {' '.join(group_members)}\n")
            f.write("set auto_zoom, on\n")
            if has_pocket:
                f.write("zoom pli_pocket_residues, 6\n")

        print(f"Protein–ligand interaction script written to {savepath_str}.")
        return savepath_str

    def _classify_ligand_features(self) -> dict:
        """
        Identify donors, acceptors, ionizable groups, aromatic rings, and hydrophobic atoms on the ligand.
        """
        donors: list[tuple[int, int, tuple[int, ...]]] = []
        acceptors: list[tuple[int, int]] = []
        hydrophobes: set[int] = set()

        for heavy_idx in range(self.ligand_no_h.GetNumAtoms()):
            full_idx = self.heavy_to_full_idx[heavy_idx]
            full_atom = self.ligand.GetAtomWithIdx(full_idx)
            heavy_atom = self.ligand_no_h.GetAtomWithIdx(heavy_idx)
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
                elif anum == 7 and chg <= 0 and full_atom.GetTotalDegree() < 4:
                    # Exclude non-accepting pyrrole-like, amide, or sulfonamide NH donors
                    has_h = bool(explicit_hs) or heavy_atom.GetTotalNumHs() > 0
                    attached_to_xo = any(
                        nb.GetAtomicNum() in (6, 16)
                        and any(
                            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(nb).GetAtomicNum() == 8
                            for b in nb.GetBonds()
                        )
                        for nb in heavy_atom.GetNeighbors()
                    )
                    if not (has_h and (heavy_atom.GetIsAromatic() or attached_to_xo)):
                        acceptors.append((heavy_idx, full_idx))
                elif anum == 16 and chg <= 0 and heavy_atom.GetDegree() <= 2:
                    acceptors.append((heavy_idx, full_idx))

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
                for match in self.ligand_no_h.GetSubstructMatches(q):
                    pos_atoms.add(match[0])

        neg_atoms: set[int] = {
            atom.GetIdx()
            for atom in self.ligand_no_h.GetAtoms()
            if atom.GetFormalCharge() < 0
        }
        if _NEG_ACID_O_QUERY is not None:
            for match in self.ligand_no_h.GetSubstructMatches(_NEG_ACID_O_QUERY):
                neg_atoms.add(match[0])

        lig_conf = self.ligand_no_h.GetConformer()
        arom_rings: list[tuple[tuple[int, ...], np.ndarray, np.ndarray]] = []
        for ring in self.ligand_no_h.GetRingInfo().AtomRings():
            if len(ring) in (5, 6) and all(
                self.ligand_no_h.GetAtomWithIdx(i).GetIsAromatic() for i in ring
            ):
                coords = np.array([lig_conf.GetAtomPosition(i) for i in ring], dtype=float)
                centroid = coords.mean(axis=0)
                normal = _ring_normal(coords)
                arom_rings.append((tuple(ring), centroid, normal))

        return {
            "donors": donors,
            "acceptors": acceptors,
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
                    prot_acceptors.append((res, res_key, p_idx, atom_name, is_bb, p_pos))

        best: dict[tuple[tuple, str, int], dict] = {}

        # Ligand Donor -> Protein Acceptor
        for lig_h_idx, lig_f_idx, _ in lig_feats["donors"]:
            d_atom = self.ligand.GetAtomWithIdx(lig_f_idx)
            d_pos = np.array(lig_conf_full.GetAtomPosition(lig_f_idx), dtype=float)
            for res, res_key, _, atom_name, is_bb, a_pos in prot_acceptors:
                dist = float(np.linalg.norm(d_pos - a_pos))
                if not (2.0 <= dist <= hbond_dist):
                    continue
                ok, ang = _eval_hbond_geometry(
                    d_atom, d_pos, a_pos, lig_conf_full, hbond_angle, DEFAULT_HBOND_HEAVY_ANGLE
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
                a_pos = np.array(lig_conf_full.GetAtomPosition(lig_f_idx), dtype=float)
                dist = float(np.linalg.norm(d_pos - a_pos))
                if not (2.0 <= dist <= hbond_dist):
                    continue
                ok, ang = _eval_hbond_geometry(
                    d_atom, d_pos, a_pos, prot_conf, hbond_angle, DEFAULT_HBOND_HEAVY_ANGLE
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
            canon = res["residue_name"]
            pdb_res = res["pdb_residue_name"]

            pos_names = set(PROTEIN_POS_ATOMS.get(canon, ()))
            if pdb_res == "HIP":
                pos_names.update({"ND1", "NE2"})

            neg_names = set(PROTEIN_NEG_ATOMS.get(canon, ()))
            if pdb_res in ("ASH", "GLH"):
                neg_names.clear()

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
                        is_parallel = dist <= 5.0 and theta <= 40.0
                        is_t_shaped = dist <= pi_stacking_dist and 50.0 <= theta <= 90.0
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
                        if dist <= cation_pi_dist:
                            key = (res_key, (l_pos_idx,))
                            if key not in best_cat or dist < best_cat[key]["distance"]:
                                best_cat[key] = _make_row(
                                    res, res_key, False, "Cation_Pi", (l_pos_idx,),
                                    "RING", dist, None, l_pos, p_cent,
                                )

            # Protein cations (LYS NZ / ARG CZ) -> Ligand aromatic rings (Cation_Pi)
            cation_names = {"NZ"} if canon == "LYS" else ({"CZ"} if canon == "ARG" else set())
            if cation_names and lig_feats["arom_rings"]:
                for p_idx in res["heavy_atom_indices"]:
                    atom_name = self.protein.GetAtomWithIdx(p_idx).GetPDBResidueInfo().GetName().strip().upper()
                    if atom_name not in cation_names:
                        continue
                    p_pos = np.array(prot_conf.GetAtomPosition(p_idx), dtype=float)
                    for l_ring, l_cent, _ in lig_feats["arom_rings"]:
                        dist = float(np.linalg.norm(p_pos - l_cent))
                        if dist <= cation_pi_dist:
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
) -> tuple[bool, Optional[float]]:
    """
    Validate hydrogen-bond angle using explicit H neighbors when present, falling back to heavy-atom
    antecedent geometry on unprotonated structures.
    """
    h_neighbors = [nbr for nbr in donor_atom.GetNeighbors() if nbr.GetAtomicNum() == 1]
    if h_neighbors:
        best_ang = -1.0
        best_ok = False
        for h_atom in h_neighbors:
            h_pos = np.array(donor_conf.GetAtomPosition(h_atom.GetIdx()), dtype=float)
            h_acc_dist = float(np.linalg.norm(h_pos - acceptor_pos))
            ang = _angle_deg(donor_pos, h_pos, acceptor_pos)
            if ang > best_ang:
                best_ang = ang
                best_ok = h_acc_dist <= 2.6 and ang >= min_h_angle
        return best_ok, best_ang

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
