import pandas as pd
import numpy as np
import py3Dmol
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdFMCS, AllChem
from rdkit.Chem.Crippen import MolLogP
import matplotlib.colors as mcolors
from typing import Optional, Union


class SAR:
    def __init__(self, data: pd.DataFrame, smi_col: str = "smiles", act_col: str = "activity", units: str = "nM"):
        self.atom_difference = None
        self.smi_col = smi_col
        self.act_col = act_col

        # check if act_col is type float
        if data[act_col].dtype == 'float':
            # print(f"Column '{act_col}' is type float")  # for debugging
            pass
        else:
            # print(f"Converting column '{act_col}' to type float")  # for debugging
            data[act_col] = data[act_col].str.replace(',', '').astype(float)

        # calculate and add pIC50 and append to data
        act_list = data[act_col].tolist()
        pIC = _calculate_pic50(act_list, units)
        data['pIC50'] = pIC
        self.data = data

    # todo add different fp types?
    def calc_LiPE(self, smi_col: Optional[list] = None, act_col: Optional[list] = None):
        """
        Calculate the Lipophilic Efficiency (LiPE) for a query molecule. The LiPE is calculated as LipE = pIC50 - LogP
        :param smi_col: Optional[list]
            Designate the smiles cols.
        :param act_col: Optional[list]
            Designate the activities cols.
        :return:
        """
        global smi_list
        if smi_col is None:
            smi_col = self.smi_col
        if act_col is None:
            act_col = self.act_col

        data = self.data
        smi_list = data[smi_col].tolist()
        pIC50_list = data['pIC50'].tolist()

        # calculate LogP
        logp = []
        for x in smi_list:
            mol = Chem.MolFromSmiles(x)
            calc_logp = MolLogP(mol)
            logp.append(calc_logp)

        # calculate LiPE
        LiPE = []
        for x, y in zip(pIC50_list, logp):
            calc_LiPE = x - y
            # result = np.round(calc_LiPE, decimals=3)
            # LiPE.append(float(result))
            LiPE.append(calc_LiPE)

        data["LiPE"] = LiPE
        self.data = data

        return self.data

    def get_sali(self, smi_col: Optional[list] = None):
        """
        Calculate the Structure-Activity Landscape Index (SALI) between a pair of molecules.
        :param smi_col: Optional[list]
            Designate the smiles cols.
        :return:
        """
        if smi_col is None:
            smi_col = self.smi_col

        data = self.data
        smi_list = data[smi_col].tolist()
        pIC50_list = data['pIC50'].tolist()

        # add fp col
        mol_list = [Chem.MolFromSmiles(x) for x in smi_list]
        data['fp'] = [Chem.RDKFingerprint(x) for x in mol_list]
        fp_list = data['fp'].tolist()

        sal_list = []
        for i, fp in enumerate(data.fp):
            sim_list = DataStructs.BulkTanimotoSimilarity(fp, fp_list)
            for j in range(0, i):
                # get mol name
                mol_name_1 = data.name.iloc[i]
                mol_name_2 = data.name.iloc[j]

                if pIC50_list[i] >= pIC50_list[j]:
                    mol_name_1, mol_name_2 = mol_name_1, mol_name_2
                else:
                    mol_name_1, mol_name_2 = mol_name_1, mol_name_2

                delta = abs(pIC50_list[i] - pIC50_list[j])
                sim = sim_list[j]
                sal_list.append(
                    [mol_name_1, mol_name_2, sim, delta / (1 - sim + 0.001), smi_list[i], smi_list[j], pIC50_list[i],
                     pIC50_list[j]])

        sal_df = pd.DataFrame(sal_list,
                              columns=["mol_1", "mol_2", "tanimoto", "SALI", "smiles_1", "smiles_2", "pIC50_1",
                                       "pIC50_2"])

        return sal_df

    def highlight_cliffs(self, smi_col: Optional[list] = None, ncols: int = 2, subsize: tuple = (400, 400),
                         legend: list = None, highlight_color: str = None, radius: int = 0.3, SVG: bool = False,
                         savepath: str = None):
        """
        Draw and highlight differing structures in 2D. For each molecule, the Maximum Common Substructure (MCS) is
        identified. Only differing functional groups will be highlighted. Information is pulled from the pd.DataFrame
        given when initializing class SAR().
        :param smi_col: str
            Column name containing SMILES strings.
        :param ncols: int
            Set the number of columns for drawn molecules.
        :param subsize: tuple
            Set the drawing size for each molecule.
        :param legend: list
            Set the legend for each moleucle in the grid.
        :param highlight_color: str
            Set the color style of the highlights. Names must be found in Matplotlib.
        :param radius: int
            Set the highlight radius.
        :param SVG: bool
            Whether to output image in SVG format. Defaults to False, giving a .png image.
        :param savepath: str
            Set the savepath for the image.
        :return:
        """
        if smi_col is None:
            smi_col = self.smi_col

        # extract smi and convert to RDKit object
        smi_list = self.data[smi_col].tolist()
        mol_list = [Chem.MolFromSmiles(smi) for smi in smi_list]

        # checks
        if len(mol_list) == 0:
            raise ValueError("No valid molecules found from provided SMILES")

            # if only one molecule, handle specially (can't find MCS)
        if len(mol_list) == 1:
            print("Only one molecule provided, no MCS can be calculated.")
            img = Draw.MolsToGridImage(
                mol_list,
                molsPerRow=1,
                subImgSize=(400, 400),
                legends=[f"Molecule {i + 1}" for i in range(len(mol_list))]
            )
            if savepath:
                img.save(savepath)
            return img

        # find mcs and create RDKit object
        mcs = rdFMCS.FindMCS(mol_list)
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

        # align mols
        try:
            ref = Chem.MolFromSmiles(smi_list[0])
            AllChem.Compute2DCoords(ref)
            for mol in mol_list:
                AllChem.GenerateDepictionMatching2DStructure(mol, ref)
        except Exception as e:
            print(f"{e}\nCannot Align Molecules.")

        # extract matching atoms
        atom_match = []
        atom_difference = []
        for mol in mol_list:
            match_atoms = mol.GetSubstructMatch(mcs_mol)
            atom_match.append(match_atoms)

            diff_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in match_atoms]
            atom_difference.append(diff_atoms)

        # save for view
        self.atom_difference = atom_difference

        if legend is None:
            try:
                legend = self.data['name'].tolist()
            except Exception as e:
                raise Exception(f"{e}\nCannot find 'name' column. No legend set!")

        if highlight_color is None:
            highlight_colors = None
        else:
            highlight = _color_to_rgb(highlight_color)
            # give color palette to each mol
            highlight_colors = [{atom_idx: highlight for atom_idx in atom_list} for atom_list in atom_difference]

        # set highlight radius
        opts = Draw.MolDrawOptions()
        opts.highlightRadius = radius

        # draw molecules
        img = Draw.MolsToGridImage(
            mols=mol_list,
            molsPerRow=ncols,
            subImgSize=subsize,
            legends=legend,
            highlightAtomLists=atom_difference,
            highlightAtomColors=highlight_colors,
            drawOptions=opts,
            useSVG=SVG,
        )

        return img

    def output_cliffs(self, mols: Chem.Mol = None, savepath: str = None):
        """
        Output the activity cliffs for molecule as .pml file for rendering in PyMOL.
        :param mols: Chem.Mol
            Input the RDKit molecule object to generate the activity cliffs.
        :param savepath: str
            Savepath for the .pml file.
        :return:
        """
        if mols is None:
            mol = self.mols
        else:
            mol = mols

        flat_indices = [idx for sublist in self.atom_difference for idx in sublist]
        # check if indices found in query mol
        flat_indices = [idx for idx in flat_indices if idx < mol.GetNumAtoms()]

        with open(savepath, "w") as f:
            # define color scheme
            f.write("set_color Hydrophob, [46, 204, 113]\n")
            f.write("set_color HDonor, [33, 150, 243]\n")
            f.write("set_color HAcceptor, [244, 67, 54]\n")
            f.write("set_color DualH, [255, 0, 255]\n")
            f.write("set_color Aromatic, [255, 235, 59]\n")

            # generate aromatic psudoatom points
            aromatic_atom_indices = set()
            ring_info = mol.GetRingInfo()
            for x, ring in enumerate(ring_info.AtomRings()):
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    if any(idx in flat_indices for idx in ring):
                        conf = mol.GetConformer()
                        coords = [conf.GetAtomPosition(idx) for idx in ring]
                        centroid = np.mean(coords, axis=0)

                        obj_name = f"Aromatic_{x}"
                        f.write(
                            f"pseudoatom {obj_name}, pos=[{centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f}]\n")
                        f.write(f"show spheres, {obj_name}\n")
                        f.write(f"hide nonbonded, {obj_name}\n")  # This removes the "plus" sign
                        f.write(f"set sphere_scale, 0.5, {obj_name}\n")
                        f.write(f"color Aromatic, {obj_name}\n")
                        f.write(f"set sphere_transparency, 0.0, {obj_name}\n")

                        for idx in ring:
                            aromatic_atom_indices.add(idx)

            # generate atom points
            for idx in flat_indices:
                atom = mol.GetAtomWithIdx(idx)
                symbol = atom.GetSymbol()
                pos = conf.GetAtomPosition(idx)

                donor = symbol in ['N', 'O', 'S'] and atom.GetTotalNumHs() > 0
                acceptor = symbol in ['N', 'O'] and atom.GetFormalCharge() <= 0

                color = None
                label = ""
                if donor and acceptor:
                    color, label = "DualH", "Dual"
                elif donor:
                    color, label = "HDonor", "Donor"
                elif acceptor:
                    color, label = "HAcceptor", "Acceptor"
                elif idx not in aromatic_atom_indices and symbol in ['C', 'Cl', 'Br', 'I']:
                    color, label = "Hydrophob", "Hydrophobic"

                if color:
                    obj_name = f"Pin_{label}_{idx}"
                    f.write(f"pseudoatom {obj_name}, pos=[{pos.x:.3f}, {pos.y:.3f}, {pos.z:.3f}]\n")
                    f.write(f"show spheres, {obj_name}\n")
                    f.write(f"color {color}, {obj_name}\n")
                    f.write(f"set sphere_transparency, 0.3, {obj_name}\n")
                    # f.write(f"hide nonbonded, {obj_name}\n")  # This removes the "plus" sign
                    # Smaller spheres for atoms inside aromatic rings
                    scale = 0.7 if idx in aromatic_atom_indices else 1.0
                    f.write(f"set sphere_scale, {scale}, {obj_name}\n")

            f.write("zoom all\n")

    def view_cliffs(self, mols: Union[Chem.Mol, list[Chem.Mol]] = None, protein_path: str = None,
                    window: tuple = (500, 500)):
        """
        View the activity cliffs of molecules in py3Dmol. Only atoms not found using Maximum Common Substructure (MCS)
        will be highlighted.
        :param mols: Union[Chem.Mol, list[Chem.Mol]]
            RDKit molecule object for rendering. Should be the same as the smiles given as the pd.DataFrame input when
            initializing SAR.
        :param protein_path: str
            Filepath to the target protein structure.
        :param window: tuple
            Set the windows size of the visualization window.
        :return:
        """

        self.mols = mols
        self.window = window
        self.protein_path = protein_path

        # dropdown menu
        if isinstance(mols, Chem.Mol):
            drop_options = [("Molecule 1", 0)]
        elif mols is None:
            raise ValueError(f"No valid RDKit molecules given!")
        else:
            drop_options = [(f"Molecule {i + 1}", i) for i in range(len(mols))]

        import ipywidgets as widgets
        dropdown = widgets.Dropdown(
            options=drop_options,
            value=0,
            description="Select:",
            style={"description_width": "initial"}
        )

        widgets.interact(self._render_cliffs, index=dropdown)

    def _render_cliffs(self, index):
        """
        Render molecules with fog effect on differing atoms
        """
        if isinstance(self.mols, Chem.Mol):
            mol = self.mols
        else:
            mol = self.mols[index]
        mol_block = Chem.MolToMolBlock(mol)

        viewer = py3Dmol.view(width=self.window[0], height=self.window[1])
        viewer.setBackgroundColor("white")
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({'stick': {'radius': 0.15, 'colorscheme': 'grayCarbon'}})
        viewer.zoomTo()

        # set protein
        if self.protein_path:
            with open(self.protein_path, 'r') as f:
                pdb_data = f.read()
            viewer.addModel(pdb_data, "pdb")
            # protein style
            viewer.setStyle({'model': -1}, {'cartoon': {'color': 'lightgray', 'opacity': 0.6},
                                            'line': {'color': 'lightgray', 'opacity': 0.3}})

        # set ligand
        viewer.addModel(Chem.MolToPDBBlock(mol), "mol")
        viewer.setStyle({'model': -1}, {'stick': {'radius': 0.15, 'colorscheme': 'grayCarbon'}})

        flat_indices = [idx for sublist in self.atom_difference for idx in sublist]
        # check if indices found in query mol
        flat_indices = [idx for idx in flat_indices if idx < mol.GetNumAtoms()]

        # track aromatic atoms
        aromatic_atom_indices = set()
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            # check if aromatic atom not in MCS
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                if any(idx in flat_indices for idx in ring):
                    # calculate centroid
                    conf = mol.GetConformer()
                    coords = [conf.GetAtomPosition(idx) for idx in ring]
                    centroid = np.mean(coords, axis=0)

                    # add sphere
                    viewer.addSphere({
                        'center': {'x': centroid[0], 'y': centroid[1], 'z': centroid[2]},
                        'radius': 0.6,
                        'color': 'gold',
                        'opacity': 1.0
                    })
                    # tag aromatic atoms
                    for idx in ring: aromatic_atom_indices.add(idx)

        # track other atoms
        for idx in flat_indices:
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()

            # tag donor/acceptor
            donor = symbol in ['N', 'O', 'S'] and atom.GetTotalNumHs() > 0
            acceptor = symbol in ['N', 'O'] and atom.GetFormalCharge() <= 0
            # use if halogens are listed as acceptors
            # acceptor = symbol in ['N', 'O', 'F', 'Cl', 'Br', 'I'] and atom.GetFormalCharge() <= 0

            # set color var
            color = None
            # donor/acceptor colorscheme
            if donor and acceptor:
                color = 'magenta'
            elif donor:
                color = 'blue'
            elif acceptor:
                color = 'red'

            # hydrophobic colorscheme
            elif idx not in aromatic_atom_indices:
                if symbol in ['C', 'F', 'Cl', 'Br', 'I'] and not atom.GetIsAromatic():
                    # only use 'C' if halogens are going to be acceptors
                    # if symbol == 'C':
                    color = '#2ecc71'
                else:
                    color = '#7f8c8d'

            # apply style
            if color:
                viewer.addStyle({'model': -1, 'index': idx}, {
                    'sphere': {
                        'color': color,
                        'opacity': 0.7,
                        'radius': 0.7 if idx in aromatic_atom_indices else 1.0
                    }
                })

        viewer.zoomTo()
        return viewer.show()


def _calculate_pic50(activity: Optional[list], units: str = "nM"):
    # check values
    for x in activity:
        if np.any(x <= 0):
            raise ValueError("Error! Input Activity must be greater than 0!")

    pic50 = []
    for x in activity:
        if units == "nM":
            converted_ic = 9 - np.log10(x)
        elif units == "uM" or units == "µM":
            converted_ic = 6 - np.log10(x)
        elif units == "mM":
            converted_ic = 3 - np.log10(x)
        else:
            raise ValueError(f"Cannot Convert {x}!")

        # converted_ic = np.round(converted_ic, decimals=3)
        pic50.append(converted_ic)

    return pic50


def _color_to_rgb(color_input):
    """Support function to convert Matplotlib colors to RGB for RDKit highlighting"""
    try:
        return mcolors.to_rgb(color_input)
    except ValueError:
        return f"Error: '{color_input}' is not a recognized color name or hex code."


if __name__ == "__main__":
    import doctest

    doctest.testmod()
