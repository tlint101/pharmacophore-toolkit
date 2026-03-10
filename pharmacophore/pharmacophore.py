import warnings
import collections
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Union, Optional

from pharmacophore.constants import feature_factory, FEATURES, FEATURE_COLORS, color_convert


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

    def read_sdf(self, sdf_file: str = None, verbose: bool = False):
        """
        Read sdf files.
        :param sdf_file: Str
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

    def to_df(self, mols: list = None, mol_name: list = None, features: Optional[Union[str, dict]] = None):
        """
        From a list of containing Chem.MOl and a list containing molecule names, create a dataframe displaying matching
        features the molecules. Method will default to "default" settings. Users can use "rdkit" to set the model to
        default RDKit features, or give a dictionary with feature:list pair, where the list contains SMARTS string for a
        given feature.
        :param mols: list
            A list containing molecules in Chem.Mol format.
        :param mol_name: list
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

    def calc_pharm(self, mol: Chem.Mol = None, features: Optional[Union[str, dict]] = None):
        """
        Generate a list of pharmacophore features and position from a molecule.
        :param mol: Chem.Mol
            A molecule in ROMol format.
        :param features: str
            Designate the type of pharmacophore features to calculate from. Will default to "default" settings. This
            param can also accept "rdkit" to use RDKit default features or give a dictionary containing feature:list
            pair, where the list contains SMARTS string for a given feature.
        :return: list
            A list of pharmacophore features, matching atom index, and position from a molecule.
        """
        global pharmacophore
        # set instance variable
        if features is None:
            features = self.features

        # calculate pharmacophore by rdkit or pharmacophore dict
        if features == 'rdkit':
            pharmacophore = self._calc_rdkit(mol)
        elif features == 'default':
            pharmacophore = self._calc_pharmacophore(mol)
        elif isinstance(features, dict):
            pharmacophore = self._calc_pharmacophore(mol, features)
        else:
            raise ValueError("Only 'default', 'rdkit', or custom features as dictionary are supported!")

        # set features to self
        self.features_list = pharmacophore

        return pharmacophore

    def add_feats(self, mol: Chem.Mol = None, substruct: str = None, type: str = None):
        """
        Optional method to add specific features if not present using default or RDKit rules.
        :param mol: Chem.Mol
            Main molecule in ROMol format.
        :param substruct: str
            A SMARTS string for specific substructure to match
        :param type: str
            Set the type of interaction for substructure. Only 'Donor', 'Acceptor', 'Aromatic', and 'Hydrophobe' is
            allowed!
        :return:
        """
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

    def output_features(self, feature_list: Optional[list] = None, savepath: str = None, type: str = "sphere",
                        sphere_size: float = 0.7, transparency: float = 0.2, color: dict = None):
        """
        Output features as a .pml format for visualization in PyMol.
        :param feature_list: Optional[list]
            A list containing features, corresponding atom number, and 3D position. Preferably genreated using the
            calc_pharm method.
        :param savepath: str = None
            Must be a file in .pml format.
        :param type: str = "Sphere"
            Set the pharmacophore representation. Can be "sphere" or "surface". If surface is used, transparency can be
            set for each pharmacophore.
        :param sphere_size: float = 0.5
            Set size of spheres.
        :param transparency: float = 0.2
            Set teh transparency of pharmacophore. Only works when the pharmacophore type is set to "surface".
        :param color: dict
            Set color for the pharmacophores. Must be given as Acceptor: color where color is a tuple for RGB or a
            string for a specific color to be translated into RGB format.
        :return:
        """
        if feature_list is None:
            feature_list = self.features_list

        with open(savepath, "w") as f:
            # feat_factory = feature_factory
            # features = feat_factory.GetFeaturesForMol(mol)
            print(f"Number of features: {len(feature_list)}")

            # set feature colors
            if color is None:
                for feat in FEATURE_COLORS:
                    f.write(f"set_color {feat}_color, {FEATURE_COLORS[feat]}\n")
            # if custom color, convert to RGB
            elif isinstance(color, dict):
                for shade in color:
                    rgb = color_convert(color[shade])
                    # print(rgb)  # for troubleshooting
                    f.write(f"set_color {shade}_color, {rgb}\n")
            else:
                warnings.warn("Issue with color features!")

            # to give sequential numbering for each group:
            feature_counts = {
                "Acceptor": 0, "Donor": 0, "Hydrophobe": 0, "Aromatic": 0, "LumpedHydrophobe": 0, "PosIonizable": 0}

            # get features
            for feat in feature_list:
                feature = feat[0]  # extract feature type
                feature_counts[feature] += 1  # give feature count
                count = feature_counts[feature]  # get current count
                # get feature position
                pos_x = feat[2]
                pos_y = feat[3]
                pos_z = feat[4]

                f.write(
                    # f"pseudoatom {feature}_{count}, pos=[{pos_x}, {pos_y}, {pos_z}], color={color_type}\n"
                    f"pseudoatom {feature}_{count}, pos=[{pos_x}, {pos_y}, {pos_z}]\n"
                )

            #todo add mol name in front of pharmacophre
            # set color and sphere size
            f.write("\n")
            if type == "sphere":
                f.write("show spheres, Acceptor_*\n")
                f.write("color acceptor_color, Acceptor_*\n")
                f.write("\n")
                f.write("show spheres, Donor_*\n")
                f.write("color donor_color, Donor_*\n")
                f.write("\n")
                f.write("show spheres, Hydrophobe_*\n")
                f.write("color hydrophobe_color, Hydrophobe_*\n")
                f.write("\n")
                f.write("show spheres, Aromatic_*\n")
                f.write("color aromatic_color, Aromatic_*\n")
                f.write("\n")
                f.write("show spheres, LumpedHydrophobe_*\n")
                f.write("color lumpedhydrophobe, LumpedHydrophobe_*\n")
                f.write("\n")
                f.write("show spheres, PosIonizable*\n")
                f.write("color posionizable, PosIonizable*\n")
                f.write("\n")
                f.write(f"set sphere_scale, {sphere_size}\n")  # Adjust sphere size in PyMOL
            elif type == "surface":
                f.write("show surface, Acceptor_*\n")
                f.write("color acceptor_color, Acceptor_*\n")
                f.write(f"set transparency, {transparency}, Acceptor_*")
                f.write("\n")
                f.write("show surface, Donor_*\n")
                f.write("color donor_color, Donor_*\n")
                f.write(f"set transparency, {transparency}, Donor_*\n")
                f.write("\n")
                f.write("show surface, Hydrophobe_*\n")
                f.write("color hydrophobe_color, Hydrophobe_*\n")
                f.write(f"set transparency, {transparency}, Hydrophobe_*\n")
                f.write("\n")
                f.write("show surface, Aromatic_*\n")
                f.write("color aromatic_color, Aromatic_*\n")
                f.write(f"set transparency, {transparency}, Aromatic_*\n")
                f.write("\n")
                f.write("show surface, LumpedHydrophobe_*\n")
                f.write("color lumpedhydrophobe, LumpedHydrophobe_*\n")
                f.write(f"set transparency, {transparency}, LumpedHydrophobe_*\n")
                f.write("\n")
                f.write("show surface, PosIonizable*\n")
                f.write("color posionizable, PosIonizable*\n")
                f.write(f"set transparency, {transparency}, PosIonizable*\n")
                f.write("\n")
            else:
                raise ValueError("Issue with type {type}! Only 'sphere' or 'surface' allowed!")

        print(f"Feature visualization script written to {savepath}.")

    def _calc_pharmacophore(self, mol: Chem.Mol = None, features: Optional[dict] = None):
        """
        Calculate pharmacophore features from a molecule using dict from constants
        :param mol: Chem.Mol
            Input molecule in ROMol format.
        :param features: Optional[dict]
            Include custom features for calculating features.
        :return:
            List of pharmacophore type and centroid coordinates.
        """
        global pharmacophore
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

        # remove duplicates for aromatic/hydrophobe set
        atom_indices = set()
        final_pharmacophore = []

        # extract aromatic/hydrophobe matches by atom index adn remove hydrophobe
        for entry in pharmacophore:
            label, atom_indexes, *rest = entry
            index = frozenset(atom_indexes)  # Convert atom indexes to a frozenset so order does not matter

            # Check if the atom indexes are already in the set
            if index not in atom_indices:
                final_pharmacophore.append(entry)
                atom_indices.add(index)
            # keep Aromatic
            elif label == 'Aromatic':
                final_pharmacophore = [e for e in final_pharmacophore if
                                       not (e[0] == 'Hydrophobe' and frozenset(e[1]) == index)]
                final_pharmacophore.append(entry)
                atom_indices.add(index)

        return final_pharmacophore

    def _calc_rdkit(self, mol: Chem.Mol = None):
        """
        Calculate pharmacophore features from a molecule using default RDKit methods.
        :param mol: Chem.Mol
            input molecule in ROMol format.
        :return:
            List of pharmacophore type and centroid coordinates.
        """
        global pharmacophore

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
def fix_bond_order(mol: Chem.Mol, template_smi: str, savepath: str = None):
    """
    Fix bond order for a given molecule. A template smiles string must be given.
    :param mol: Chem.Mol
        Molecule to fix. Must be in ROMol format.
    :param template_smi: str
        A smiles string for molecule to fix. Recommend canonical smiles string if possible.
    :param savepath: str
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


def find_matches(mol: Chem.Mol = None, patterns: list[Chem.Mol] = None, verbose=True):
    """
    Support function to visualize matches between query molecule and features.
    :param mol: Chem.Mol
        Query molecule. Must be in ROMol format.
    :param patterns: list[Chem.Mol]
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
