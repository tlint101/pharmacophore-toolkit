import math
import os
import warnings
import collections
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Union, Optional

from pharmacophore.constants import (
    feature_factory,
    FEATURES,
    FEATURE_COLORS,
    color_convert,
    format_pymol_split_spheres_block,
    resolve_pml_output_paths,
    sanitize_pymol_name,
)


def _extract_mol_name(mol: Chem.Mol, idx: int = 0) -> str:
    """Return a molecule's ``_Name`` property when non-empty, else ``Mol_{idx}``."""
    prop_name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    return prop_name or f"Mol_{idx}"


def _is_single_feature_entry(item) -> bool:
    """Return True if ``item`` matches ``[family, atom_indices, x, y, z]``."""
    return isinstance(item, (list, tuple)) and len(item) == 5 and isinstance(item[0], str)


class Pharmacophore:
    def __init__(self, features: Union[str, dict] = 'default'):
        """
        Initialize the Pharmacophore class.
        :param features: Union[str, dict]
            Set features to use. Will default to 'default'. Can also accept RDKit features or a dictionary of features.
        """
        self.sdf = None
        self.features = features
        self.features_list = None
        self.mol_names: Optional[list[str]] = None
        self._features_is_batch: bool = False

    def read_sdf(self, sdf_file: Optional[str] = None, verbose: bool = False):
        """
        Read sdf files.
        :param sdf_file: Optional[str]
            File path to .sdf file.
        :param verbose: Bool
            Output description when reading .sdf file.
        :return: list
            A list of molecules in Chem.Mol format.
        """
        supplier = Chem.SDMolSupplier(sdf_file)
        # extract ROMols into a list
        mol_list = []
        if verbose:
            for mol in tqdm(supplier, desc="Reading Molecules"):
                mol_list.append(mol)
        else:
            for mol in supplier:
                mol_list.append(mol)

        self.sdf = mol_list

        return mol_list

    def feature_types(self, features: Optional[Union[str, dict]] = None):
        """
        Output feature types that will be displayed in the pharmacophore model. This will default to 'default' settings.
        However, default RDKit features can also be used or a dictionary containing a dictionary of features with SMARTS
        string can be used.
        :param features: Optional[Union[str, dict]]
            The feature type to be used in the pharmacophore model. This will default to "default" settings. Using
            "rdkit" will set the model to default RDKit features. Custom features can be used by giving a dictionary
            with feature:list pair, where the list contains SMARTS string for a given feature.

        :return:
        """
        # set instance variable
        if features is None:
            features = self.features

        global phrase
        if features == "default":
            phrase = f"Default features: \n('Donor', 'Acceptor', 'Aromatic', 'Hydrophobe')"
        elif features == "rdkit":
            phrase = f"Default features from RDKit: \n{feature_factory.GetFeatureFamilies()}"
        elif isinstance(features, dict):
            phrase = f"Custom features: \n{list(features.keys())}"

        return phrase

    def to_df(self, mols: Optional[list] = None, mol_name: Optional[list] = None,
              features: Optional[Union[str, dict]] = None):
        """
        From a list of containing Chem.MOl and a list containing molecule names, create a dataframe displaying matching
        features the molecules. Method will default to "default" settings. Users can use "rdkit" to set the model to
        default RDKit features, or give a dictionary with feature:list pair, where the list contains SMARTS string for a
        given feature.
        :param mols: Optional[list]
            A list containing molecules in Chem.Mol format.
        :param mol_name: Optional[list]
            A list containing names of molecules.
        :param features: Union[str, dict]
            The feature type to be used in the pharmacophore model. This will default to "default" settings. Using
            "rdkit" will set the model to default RDKit features. Custom features can be used by giving a dictionary
            with feature:list pair, where the list contains SMARTS string for a given feature.
        :return:
        """
        global feat_factory
        # set instance variable
        if features is None:
            features = self.features

        molecule_feature_frequencies = []

        # use default features
        if features == 'rdkit':
            feat_factory = feature_factory
            # get default features, count frequency, and append into list
            for mol in mols:
                feats = [feature.GetFamily() for feature in feat_factory.GetFeaturesForMol(mol)]
                feature_frequency = collections.Counter(feats)
                molecule_feature_frequencies.append(feature_frequency)

        # use default features
        elif features == 'default':
            feature_list = []
            # get features from constant file, count pharmacophore type, append to list
            for mol in mols:
                feats = self._calc_pharmacophore(mol)
                for feat in feats:
                    feature_list.append(feat[0])
                feature_frequency = collections.Counter(feature_list)
                molecule_feature_frequencies.append(feature_frequency)
                feature_list = []  # reset list to avoid double counting

        # use custom features
        elif isinstance(features, dict):
            # reset instance variable
            self.features = features

            feature_list = []
            # get features from constant file, count pharmacophore type, append to list
            for mol in mols:
                feats = self._calc_pharmacophore(mol, features)
                for feat in feats:
                    feature_list.append(feat[0])
                feature_frequency = collections.Counter(feature_list)
                molecule_feature_frequencies.append(feature_frequency)
                feature_list = []  # reset list to avoid double counting

        elif isinstance(features, dict) is False:
            raise ValueError("Only 'default', 'rdkit', or custom features as dictionary are supported!")

        # check if mol_name is given
        if mol_name is None:
            mol_name = []
            for x in range(len(mols)):
                name = "Mol" + str(x)
                mol_name.append(name)

        # # for troubleshooting
        # print(molecule_feature_frequencies)

        # reformat table and rename columns to molecule list
        df = pd.DataFrame(molecule_feature_frequencies)
        feature_frequencies_df = df.transpose()
        feature_frequencies_df = feature_frequencies_df.fillna(0).astype(int)
        feature_frequencies_df = feature_frequencies_df.set_axis(mol_name, axis=1)

        return feature_frequencies_df

    def _calc_single_pharm(self, mol: Chem.Mol, features: Union[str, dict]) -> list:
        """Compute pharmacophore features for a single RDKit molecule."""
        if features == 'rdkit':
            return self._calc_rdkit(mol)
        if features == 'default':
            return self._calc_pharmacophore(mol)
        if isinstance(features, dict):
            return self._calc_pharmacophore(mol, features)
        raise ValueError("Only 'default', 'rdkit', or custom features as dictionary are supported!")

    def calc_pharm(
        self,
        mol: Optional[Union[Chem.Mol, list[Chem.Mol], tuple[Chem.Mol, ...]]] = None,
        features: Optional[Union[str, dict]] = None,
    ):
        """
        Generate pharmacophore features and 3D positions from one or more molecules.
        :param mol: Optional[Union[Chem.Mol, list[Chem.Mol], tuple[Chem.Mol, ...]]]
            A single molecule or a sequence of molecules in RDKit ``Chem.Mol`` format.
        :param features: Optional[Union[str, dict]]
            Designate the type of pharmacophore features to calculate from. Will default to "default" settings. This
            param can also accept "rdkit" to use RDKit default features or give a dictionary containing feature:list
            pair, where the list contains SMARTS string for a given feature.
        :return: list
            For a single molecule, a list of ``[family, atom_indices, x, y, z]`` entries. For a sequence of
            molecules, a list of per-molecule feature lists matching the input order.
        """
        if mol is None:
            raise ValueError("mol cannot be None.")
        if features is None:
            features = self.features

        if isinstance(mol, Chem.Mol):
            pharmacophore = self._calc_single_pharm(mol, features)
            self.features_list = pharmacophore
            self.mol_names = [_extract_mol_name(mol, 0)]
            self._features_is_batch = False
            return pharmacophore

        if isinstance(mol, (list, tuple)):
            if not mol:
                raise ValueError("mol sequence cannot be empty.")
            for idx, m in enumerate(mol):
                if not isinstance(m, Chem.Mol):
                    raise TypeError(f"Expected RDKit Chem.Mol at mol[{idx}], got {type(m).__name__}.")
            all_pharmacophores = [self._calc_single_pharm(m, features) for m in mol]
            self.features_list = all_pharmacophores
            self.mol_names = [_extract_mol_name(m, idx) for idx, m in enumerate(mol)]
            self._features_is_batch = True
            return all_pharmacophores

        raise TypeError("mol must be an RDKit Chem.Mol or a sequence of Chem.Mol objects.")

    def add_feats(self, mol: Optional[Chem.Mol] = None, substruct: Optional[str] = None, type: Optional[str] = None):
        """
        Optional method to add specific features if not present using default or RDKit rules.
        :param mol: Optional[Chem.Mol]
            Main molecule in ROMol format.
        :param substruct: Optional[str]
            A SMARTS string for specific substructure to match
        :param type: Optional[str]
            Set the type of interaction for substructure. Only 'Donor', 'Acceptor', 'Aromatic', and 'Hydrophobe' is
            allowed!
        :return:
        """
        if self._features_is_batch:
            raise ValueError(
                "add_feats() operates on a single molecule's feature list. Call calc_pharm() on a single molecule first."
            )
        if self.features_list is None:
            self.features_list = []

        # convert substruct smarts string into ROMol
        query = [Chem.MolFromSmarts(substruct)]

        # find matches
        matches = find_matches(mol, query)

        # check type
        if type == 'Donor' or type == 'Acceptor' or type == 'Aromatic' or type == 'Hydrophobe':
            pass
        else:
            raise ValueError(f"Type {type} not supported! Only 'Donor', 'Acceptor', 'Aromatic', and 'Hydrophobe' is "
                             f"allowed!")

        # obtain specific data from each match
        for match in matches:
            result = [type, match[0], match[1][0], match[1][1], match[1][2]]
            # add custom feats to self.features_list
            self.features_list.append(result)

        return self.features_list

    def _normalize_feature_batches(
        self,
        feature_list: Optional[list],
        mol_names: Optional[Union[str, list[str], tuple[str, ...]]],
    ) -> tuple[list[list], list[str], bool]:
        """Normalize single or multi-molecule feature inputs and resolve per-molecule names."""
        if feature_list is None:
            if self.features_list is None:
                raise ValueError("No pharmacophore features available. Call calc_pharm() or pass feature_list.")
            if self._features_is_batch:
                batches = [list(feats) for feats in self.features_list]
            else:
                batches = [list(self.features_list)]
            if self.mol_names is not None and len(self.mol_names) == len(batches):
                resolved_names = list(self.mol_names)
            else:
                resolved_names = [f"Mol_{i}" for i in range(len(batches))]
        else:
            if not isinstance(feature_list, (list, tuple)):
                raise ValueError("feature_list must be a list of features or a list of per-molecule feature lists.")
            if len(feature_list) == 0:
                batches = [[]]
            elif all(_is_single_feature_entry(x) for x in feature_list):
                batches = [list(feature_list)]
            elif all(
                isinstance(x, (list, tuple)) and all(_is_single_feature_entry(entry) for entry in x)
                for x in feature_list
            ):
                batches = [list(x) for x in feature_list]
            else:
                raise ValueError(
                    "Invalid feature_list format. Expected [family, atom_indices, x, y, z] entries "
                    "or a sequence of per-molecule feature lists."
                )

            if (
                feature_list is self.features_list
                and self.mol_names is not None
                and len(self.mol_names) == len(batches)
            ):
                resolved_names = list(self.mol_names)
            else:
                resolved_names = [f"Mol_{i}" for i in range(len(batches))]

        if mol_names is not None:
            if isinstance(mol_names, str):
                if len(batches) != 1:
                    raise ValueError(
                        f"Expected {len(batches)} molecule names in mol_names, got single string {mol_names!r}."
                    )
                resolved_names = [mol_names]
            elif isinstance(mol_names, (list, tuple)):
                if len(mol_names) != len(batches):
                    raise ValueError(
                        f"Length of mol_names ({len(mol_names)}) must match number of molecules ({len(batches)})."
                    )
                resolved_names = ["" if n is None else str(n) for n in mol_names]
            else:
                raise TypeError("mol_names must be a string or a sequence of strings.")

        return batches, resolved_names, len(batches) > 1

    def output_features(
        self,
        feature_list: Optional[list] = None,
        mol_names: Optional[Union[str, list[str]]] = None,
        type: str = "sphere",
        sphere_size: float = 0.7,
        transparency: float = 0.2,
        color: Optional[dict] = None,
        savepath: Optional[Union[str, os.PathLike]] = None,

    ) -> Union[str, list[str]]:
        """
        Output features for one or more molecules as ``.pml`` script(s) for visualization in PyMOL.

        For a single molecule, ``savepath`` may be a ``.pml`` file path or an output directory
        (writing ``<dir>/<mol_name>_pharmacophore.pml``). For multiple molecules, ``savepath`` must be an
        output directory and writes ``<mol_name>_pharmacophore.pml`` per molecule with molecule-prefixed
        PyMOL object and group names.
        :param feature_list: Optional[list]
            A feature list from ``calc_pharm(mol)``, a list of per-molecule feature lists from
            ``calc_pharm(mols)``, or ``None`` to use the last result stored by ``calc_pharm()``.
        :param mol_names: Optional[Union[str, list[str]]] = None
            Optional molecule name(s) used for ``<mol_name>_pharmacophore.pml`` filenames and PyMOL object prefixes.
            Defaults to names recorded by ``calc_pharm()`` or ``Mol_0``, ``Mol_1``, etc.
        :param type: str = "sphere"
            Set the pharmacophore representation. Can be "sphere" or "surface". If surface is used, transparency can be
            set for each pharmacophore.
        :param sphere_size: float = 0.7
            Set size of spheres.
        :param transparency: float = 0.2
            Set the transparency of pharmacophore. Only works when the pharmacophore type is set to "surface".
        :param color: Optional[dict]
            Set color for the pharmacophores. Must be given as ``{family: color}`` where color is a tuple for RGB
            or a string for a specific color to be translated into RGB format.
        :param savepath: Optional[Union[str, os.PathLike]] = None
            Required ``.pml`` file path (single molecule) or output directory path (single or multiple molecules).
        :return: Union[str, list[str]]
            Written ``.pml`` path for a single molecule, or list of written paths for multiple molecules.
        """
        if type not in ("sphere", "surface"):
            raise ValueError(f"Issue with type {type!r}! Only 'sphere' or 'surface' allowed!")

        batches, resolved_names, is_multi = self._normalize_feature_batches(feature_list, mol_names)
        out_paths, safe_names = resolve_pml_output_paths(savepath, resolved_names, suffix="pharmacophore")

        for feats, out_path, mol_prefix in zip(batches, out_paths, safe_names):
            _write_single_features_pml(
                feature_list=feats,
                mol_prefix=mol_prefix,
                savepath=out_path,
                type=type,
                sphere_size=sphere_size,
                transparency=transparency,
                color=color,
                verbose=not is_multi,
            )

        if is_multi:
            print(f"Number of molecules: {len(out_paths)}")
            print(f"Feature visualization scripts written to {str(savepath).strip()}.")

        return out_paths if is_multi else out_paths[0]

    def _calc_pharmacophore(self, mol: Optional[Chem.Mol] = None, features: Optional[dict] = None):
        """
        Calculate pharmacophore features from a molecule using dict from constants
        :param mol: Optional[Chem.Mol]
            Input molecule in ROMol format.
        :param features: Optional[dict]
            Include custom features for calculating features.
        :return:
            List of pharmacophore type and centroid coordinates.
        """
        # read in custom feature dict
        if features is None:
            constant_feats = FEATURES
        else:
            constant_feats = features

        # hold matches
        matches = {}

        # Identify matches to feature dict
        for key, value in constant_feats.items():
            try:
                query = [Chem.MolFromSmarts(smarts) for smarts in value]
                matches[key] = find_matches(mol, query, verbose=False)
            except:
                pass
        # remove duplicate SMARTS matches
        cleaned_matches = {}
        for key, value in matches.items():
            unique_lists = []
            for list in value:
                if list not in unique_lists:
                    unique_lists.append(list)
            cleaned_matches[key] = unique_lists
        # create list containing pharmacophore and centroid coordinates
        pharmacophore = []
        for key, value in cleaned_matches.items():
            for match in value:
                # extarct feature type, atom match, and XYZ position
                p = [key, match[0], match[1][0], match[1][1], match[1][2]]
                pharmacophore.append(p)

        # Deduplicate exact (family, atom_indices) matches (preserving dual Donor + Acceptor on the same atom)
        seen_family_indices = set()
        deduped = []
        for entry in pharmacophore:
            label, atom_indexes, *rest = entry
            key = (label, frozenset(atom_indexes))
            if key not in seen_family_indices:
                seen_family_indices.add(key)
                deduped.append(entry)

        # Prefer Aromatic over Hydrophobe when both match the exact same ring
        aromatic_rings = {frozenset(e[1]) for e in deduped if e[0] == "Aromatic"}
        final_pharmacophore = [
            e for e in deduped
            if not (e[0] == "Hydrophobe" and frozenset(e[1]) in aromatic_rings)
        ]

        return final_pharmacophore

    def _calc_rdkit(self, mol: Optional[Chem.Mol] = None):
        """
        Calculate pharmacophore features from a molecule using default RDKit methods.
        :param mol: Optional[Chem.Mol]
            input molecule in ROMol format.
        :return:
            List of pharmacophore type and centroid coordinates.
        """
        # load default RDKit feature factory file
        feat_factory = feature_factory

        # get list of features for molecule
        features = feat_factory.GetFeaturesForMol(mol)

        # store pharmacophore features as list
        pharmacophore = []
        for feature in features:
            fam = feature.GetFamily()
            pos = feature.GetPos()
            atom_indices = feature.GetAtomIds()
            pharmacophore_item = [fam, atom_indices, pos[0], pos[1], pos[2]]
            pharmacophore.append(pharmacophore_item)

        return pharmacophore


# Support function to fix bond order of molecule
def _write_single_features_pml(
    feature_list: list,
    mol_prefix: str,
    savepath: str,
    type: str = "sphere",
    sphere_size: float = 0.7,
    transparency: float = 0.2,
    color: Optional[dict] = None,
    verbose: bool = True,
) -> None:
    """Write a PyMOL ``.pml`` script for one molecule's pharmacophore features."""
    with open(savepath, "w") as f:
        if verbose:
            print(f"Number of features: {len(feature_list)}")
        f.write("# Drag-and-drop PyMOL script generated by pharmacophore-toolkit (Pharmacophore)\n")
        f.write("set auto_zoom, off\n\n")

        resolved_colors = {k: tuple(float(c) for c in v[:3]) for k, v in FEATURE_COLORS.items()}
        if color is not None:
            if isinstance(color, dict):
                for shade, val in color.items():
                    rgb = (
                        tuple(float(c) for c in val[:3])
                        if isinstance(val, (tuple, list)) and len(val) in (3, 4)
                        else tuple(color_convert(val))
                    )
                    resolved_colors[shade] = rgb
            else:
                warnings.warn("Issue with color features!")

        for feat in feature_list:
            fam = str(feat[0])
            if fam not in resolved_colors:
                resolved_colors[fam] = (0.5, 0.5, 0.5)

        emitted_color_names: set[str] = set()
        family_color_name: dict[str, str] = {}
        for fam, rgb in resolved_colors.items():
            safe_fam = sanitize_pymol_name(fam, fallback="Feature")
            cname = f"{mol_prefix}_{safe_fam.lower()}_color"
            family_color_name[fam] = cname
            if cname not in emitted_color_names:
                emitted_color_names.add(cname)
                f.write(f"set_color {cname}, [{rgb[0]:.4f}, {rgb[1]:.4f}, {rgb[2]:.4f}]\n")

        feature_counts: dict[str, int] = {}
        coord_groups: dict[tuple[float, float, float], list[int]] = {}
        for idx, feat in enumerate(feature_list):
            coord_key = (round(float(feat[2]), 3), round(float(feat[3]), 3), round(float(feat[4]), 3))
            coord_groups.setdefault(coord_key, []).append(idx)

        single_families: list[tuple[str, str]] = []
        seen_single_families: set[str] = set()
        group_members: list[str] = []
        split_specs = []

        for idx, feat in enumerate(feature_list):
            feature = str(feat[0])
            safe_feat = sanitize_pymol_name(feature, fallback="Feature")
            feature_counts[safe_feat] = feature_counts.get(safe_feat, 0) + 1
            count = feature_counts[safe_feat]
            pos_x, pos_y, pos_z = feat[2], feat[3], feat[4]
            obj_name = f"{mol_prefix}_{safe_feat}_{count}"
            group_members.append(obj_name)
            coord_key = (round(float(pos_x), 3), round(float(pos_y), 3), round(float(pos_z), 3))
            group = coord_groups[coord_key]
            if len(group) == 1:
                f.write(f"pseudoatom {obj_name}, pos=[{pos_x}, {pos_y}, {pos_z}]\n")
                if safe_feat not in seen_single_families:
                    seen_single_families.add(safe_feat)
                    single_families.append((safe_feat, family_color_name[feature]))
            else:
                rank = group.index(idx)
                k = len(group)
                phi0 = 2.0 * math.pi * rank / k
                phi1 = 2.0 * math.pi * (rank + 1) / k
                rgb = resolved_colors.get(feature, (0.5, 0.5, 0.5))
                split_specs.append((obj_name, float(pos_x), float(pos_y), float(pos_z), rgb, phi0, phi1))

        f.write("\n")
        for safe_feat, cname in single_families:
            pattern = f"{mol_prefix}_{safe_feat}_*"
            if type == "sphere":
                f.write(f"show spheres, {pattern}\n")
                f.write(f"color {cname}, {pattern}\n")
                f.write(f"set sphere_scale, {sphere_size}, {pattern}\n\n")
            elif type == "surface":
                f.write(f"show surface, {pattern}\n")
                f.write(f"color {cname}, {pattern}\n")
                f.write(f"set transparency, {transparency}, {pattern}\n\n")

        if split_specs:
            f.write(
                format_pymol_split_spheres_block(
                    split_specs,
                    sphere_scale=(sphere_size if type == "sphere" else None),
                    transparency=(0.0 if type == "sphere" else transparency),
                )
            )

        if group_members:
            f.write(f"group {mol_prefix}_pharmacophore, {' '.join(group_members)}\n")
        f.write("set auto_zoom, on\n")

    if verbose:
        print(f"Feature visualization script written to {savepath}.")


def fix_bond_order(mol: Chem.Mol, template_smi: str, savepath: Optional[str] = None):
    """
    Fix bond order for a given molecule. A template smiles string must be given.
    :param mol: Chem.Mol
        Molecule to fix. Must be in ROMol format.
    :param template_smi: str
        A smiles string for molecule to fix. Recommend canonical smiles string if possible.
    :param savepath: Optional[str]
        Set the save path for the molecule. Output will be in sdf format.
    :return:
    """
    # remove hydrogen if present
    mol_prep = Chem.RemoveHs(mol)

    # read template
    template = Chem.MolFromSmiles(template_smi)

    # assign bond orders to crystal ligand
    fixed_mol = AllChem.AssignBondOrdersFromTemplate(mol_prep, template)

    # Map atoms from the template to the docked pose
    template_mapping = mol_prep.GetSubstructMatch(template)

    if not template_mapping:
        raise ValueError("Failed to match template to docked pose.")

    # Transfer 3D coordinates from docked_pose to fixed_mol
    new_conformer = Chem.Conformer(fixed_mol.GetNumAtoms())
    original_conformer = mol_prep.GetConformer()

    # set coordinates of heavy atoms to fixed_mol
    for new_idx, origin_idx in enumerate(template_mapping):
        pos = original_conformer.GetAtomPosition(origin_idx)
        new_conformer.SetAtomPosition(new_idx, pos)

    # Add 3D coordinates to the molecule
    fixed_mol.AddConformer(new_conformer, assignId=True)

    # save file
    Chem.MolToMolFile(fixed_mol, savepath)


def find_matches(mol: Optional[Chem.Mol] = None, patterns: Optional[list[Chem.Mol]] = None, verbose=True):
    """
    Support function to visualize matches between query molecule and features.
    :param mol: Optional[Chem.Mol]
        Query molecule. Must be in ROMol format.
    :param patterns: Optional[list[Chem.Mol]]
        A list of molecules to match against the query molecule.
    :param verbose: Bool
        Output messages for matches.
    :return:
    """
    matches = []
    for pattern in patterns:
        matched = mol.GetSubstructMatches(pattern)
        # get centroid and coordinates for each match
        for m in matched:
            try:
                centroid = _compute_match_centroid(mol, m)
                matches.append([m, centroid])
            except:
                pass
        # output statement if no matches
        if verbose is True:
            if len(matches) == 0:
                output_message = Chem.MolToSmarts(pattern)
                print(f"No Matches to {output_message}!")
                return matches

    return matches


def _compute_match_centroid(mol, matched_pattern):
    """
    Support function to calculate centroid of matches between query molecule and features.
    :param mol:
    :param matched_pattern:
    :return:
    """
    conf = mol.GetConformer()
    positions = [conf.GetAtomPosition(i) for i in matched_pattern]
    center = np.mean(positions, axis=0)
    # convert result to float
    center = center.tolist()
    return tuple(center)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
