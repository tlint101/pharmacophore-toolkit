import json
import pandas as pd
import numpy as np
import py3Dmol
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdFMCS, AllChem, rdFingerprintGenerator, MACCSkeys
from rdkit.Chem.Crippen import MolLogP
import matplotlib.colors as mcolors
from typing import Optional, Union


class SAR:
    def __init__(self, data: pd.DataFrame, smi_col: str = "smiles", act_col: str = "activity", units: str = "nM",
                 type: str = "jupyter"):
        """
        Shared functions for the SAR() class.
        :param data: pd.DataFrame
            A DataFrame containing the SMILES and activity of each molecule.
        :param smi_col: str
            Column name containing SMILES strings.
        :param act_col: str
            Column name containing the activity values.
        :param units: str
            Units of the activity values. Accepts 'nM', 'uM' or 'mM'.
        :param type: str
            Set the output for Jupyter or Marimo notebooks. Defaults to Jupyter.
        """
        # check type
        if type.lower() not in ("jupyter", "marimo"):
            raise ValueError("Only 'jupyter' or 'marimo' accepted!")
        self.type = type.lower()

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

    def get_sali(self, smi_col: Optional[list] = None, bits: int = 1024, radius: int = 2, type: str = 'morgan'):
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

        # todo fix

        if type == 'morgan':
            print(Warning("Currently only RDKFingerprint will be used by default"))

        # add fp col
        mol_list = [Chem.MolFromSmiles(x) for x in smi_list]
        data['fp'] = [Chem.RDKFingerprint(x) for x in
                      mol_list]  # todo a function to use different types of fingerprints
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

    def highlight_cliffs(self, smi_list: Optional[list] = None, ncols: int = 2, subsize: tuple = (400, 400),
                         legend: Optional[list] = None, highlight_color: Optional[str] = None, radius: int = 0.3,
                         SVG: bool = False, savepath: Optional[str] = None):
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
        :param legend: Optional[list]
            Set the legend for each moleucle in the grid.
        :param highlight_color: Optional[str]
            Set the color style of the highlights. Names must be found in Matplotlib.
        :param radius: int
            Set the highlight radius.
        :param SVG: bool
            Whether to output image in SVG format. Defaults to False, giving a .png image.
        :param savepath: Optional[str]
            Set the savepath for the image.
        :return:
        """
        # extract smi and convert to RDKit object
        if smi_list is None:
            smi_col = self.smi_col
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

    def output_cliffs(self, mols: Optional[Chem.Mol] = None, savepath: Optional[str] = None):
        """
        Output the activity cliffs for molecule as .pml file for rendering in PyMOL.
        :param mols: Optional[Chem.Mol]
            Input the RDKit molecule object to generate the activity cliffs.
        :param savepath: Optional[str]
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

    def view_cliffs(self, mols: Union[Chem.Mol, list[Chem.Mol]] = None, protein_path: Optional[str] = None,
                    window: tuple = (500, 500), prefix: str = "cliffs"):
        """
        View the activity cliffs of molecules in py3Dmol, in either a Jupyter or Marimo notebook. Only atoms not found
        using Maximum Common Substructure (MCS) will be highlighted.
        :param mols: Union[Chem.Mol, list[Chem.Mol]]
            RDKit molecule object for rendering. Should be the same as the smiles given as the pd.DataFrame input when
            initializing SAR.
        :param protein_path: Optional[str]
            Filepath to the target protein structure.
        :param window: tuple
            Set the windows size of the visualization window.
        :param prefix: str
            Set the prefix for the saved image. Only works for Marimo notebooks.
        :return:
        """
        if mols is None:
            raise ValueError(f"No valid RDKit molecules given!")

        self.mols = [mols] if isinstance(mols, Chem.Mol) else mols
        self.window = window
        self.protein_path = protein_path

        if self.type == "jupyter":
            return self._jupyter_cliffs()
        elif self.type == "marimo":
            return self._marimo_cliffs(prefix)
        else:
            raise ValueError("Only 'jupyter' or 'marimo' accepted!")

    def _jupyter_cliffs(self):
        """Output Jupyter interactive window."""
        import ipywidgets as widgets
        dropdown = widgets.Dropdown(
            options=[(f"Molecule {i + 1}", i) for i in range(len(self.mols))],
            value=0,
            description="Select:",
            style={"description_width": "initial"}
        )

        widgets.interact(self._render_cliffs, index=dropdown)

    def _cliff_shapes(self, mol):
        """
        Support function for _render_cliffs and _marimo_cliffs. Computes the highlight spheres and the per-atom sphere
        styles for the atoms of one molecule that fall outside the MCS.
        """
        flat_indices = [idx for sublist in self.atom_difference for idx in sublist]
        # check if indices found in query mol
        flat_indices = [idx for idx in flat_indices if idx < mol.GetNumAtoms()]

        # track aromatic atoms
        spheres = []
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
                    spheres.append({
                        'center': {'x': float(centroid[0]), 'y': float(centroid[1]), 'z': float(centroid[2])},
                        'radius': 0.6,
                        'color': 'gold',
                        'opacity': 1.0
                    })
                    # tag aromatic atoms
                    for idx in ring: aromatic_atom_indices.add(idx)

        # track other atoms
        atoms = []
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

            # collect style
            if color:
                atoms.append({
                    'index': idx,
                    'color': color,
                    'opacity': 0.7,
                    'radius': 0.7 if idx in aromatic_atom_indices else 1.0
                })

        return spheres, atoms

    def _render_cliffs(self, index):
        """
        Render molecules with fog effect on differing atoms
        """
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

        spheres, atoms = self._cliff_shapes(mol)
        for sphere in spheres:
            viewer.addSphere(sphere)

        # apply style
        for atom in atoms:
            viewer.addStyle({'model': -1, 'index': atom['index']}, {
                'sphere': {'color': atom['color'], 'opacity': atom['opacity'], 'radius': atom['radius']}
            })

        viewer.show()

    def _marimo_cliffs(self, prefix):
        """
        Output Marimo interactive window.
        """
        try:
            import marimo as mo
        except Exception as e:
            return e

        # serialize molecule data
        data = []
        for i, mol in enumerate(self.mols):
            spheres, atoms = self._cliff_shapes(mol)
            data.append({
                "name": f"Molecule {i + 1}",
                "mol": Chem.MolToMolBlock(mol),
                "pdb": Chem.MolToPDBBlock(mol),
                "spheres": spheres,
                "atoms": atoms,
            })

        protein = None
        if self.protein_path:
            with open(self.protein_path, 'r') as f:
                protein = f.read()

        width, height = self.window

        html = f"""
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <div style="font:13px sans-serif;display:flex;gap:8px;align-items:center;margin-bottom:8px">
          <select id="sel"></select>
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
        const DATA = {json.dumps(data)};
        const PROTEIN = {json.dumps(protein)};
        const box = document.getElementById("viewer");
        const sel = document.getElementById("sel");
        DATA.forEach((d, i) => sel.add(new Option(d.name, i)));
        const viewer = $3Dmol.createViewer(box, {{backgroundColor: "white"}});

        function render(i) {{
          const d = DATA[i];
          viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllLabels();
          viewer.addModel(d.mol, "mol");
          viewer.setStyle({{stick: {{radius: 0.15, colorscheme: "grayCarbon"}}}});
          viewer.zoomTo();
          if (PROTEIN) {{
            viewer.addModel(PROTEIN, "pdb");
            viewer.setStyle({{model: -1}}, {{cartoon: {{color: "lightgray", opacity: 0.6}},
                                            line: {{color: "lightgray", opacity: 0.3}}}});
          }}
          viewer.addModel(d.pdb, "mol");
          viewer.setStyle({{model: -1}}, {{stick: {{radius: 0.15, colorscheme: "grayCarbon"}}}});
          d.spheres.forEach(s => viewer.addSphere(s));
          d.atoms.forEach(a => viewer.addStyle({{model: -1, index: a.index}},
            {{sphere: {{color: a.color, opacity: a.opacity, radius: a.radius}}}}));
          viewer.render();
        }}

        function snapshot() {{
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
          a.download = "{prefix}_" + DATA[+sel.value].name.replace(/\\s+/g, "_") + ".png";
          document.body.appendChild(a); a.click(); a.remove();
        }}

        document.querySelectorAll("[data-pan]").forEach(b => b.onclick = () => {{
          const [dx, dy] = b.dataset.pan.split(",").map(Number);
          viewer.translate(dx, dy, 200);
        }});
        document.getElementById("recenter").onclick = () => viewer.zoomTo();
        sel.onchange = e => render(+e.target.value);
        document.getElementById("shot").onclick = snapshot;
        render(0);
        </script>"""

        return mo.iframe(html, height=f"{height + 110}px")


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


# todo add additional fingerprints
def _calculate_fingerprint(mol, fp_type='morgan', radius=2, fpSize=1024):
    """
    Calculates molecular fingerprints for various RDKit algorithms.

    Parameters:
        mol: RDKit Mol object
        fp_type: Type of fingerprint ('morgan', 'rdkit', 'atompair', 'torsion', 'maccs')
        radius: Radius for Morgan (ignored by others)
        fpSize: Size of the bit vector (ignored by MACCS)
    """
    # 1. Map string names to their respective Generator factory functions
    generators = {
        'morgan': lambda: rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fpSize),
        'rdkit': lambda: rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=fpSize),
        'atompair': lambda: rdFingerprintGenerator.GetAtomPairGenerator(fpSize=fpSize),
        'torsion': lambda: rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=fpSize),
    }

    # 2. Handle MACCS Keys separately (they have a fixed size of 167 bits)
    if fp_type.lower() == 'maccs':
        return MACCSkeys.GenMACCSKeys(mol)

    # 3. Get the generator and produce the fingerprint
    if fp_type.lower() in generators:
        gen = generators[fp_type.lower()]()
        return gen.GetFingerprint(mol)
    else:
        raise ValueError(f"Unknown fingerprint type: {fp_type}. "
                         f"Choose from {list(generators.keys()) + ['maccs']}")


def _color_to_rgb(color_input):
    """Support function to convert Matplotlib colors to RGB for RDKit highlighting"""
    try:
        return mcolors.to_rgb(color_input)
    except ValueError:
        return f"Error: '{color_input}' is not a recognized color name or hex code."


if __name__ == "__main__":
    import doctest

    doctest.testmod()
