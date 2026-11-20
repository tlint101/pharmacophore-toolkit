import json
import os
import pandas as pd
import numpy as np
import py3Dmol
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdFMCS, AllChem, rdFingerprintGenerator, MACCSkeys
from rdkit.Chem.Crippen import MolLogP
from rdkit.Avalon import pyAvalonTools
import matplotlib.colors as mcolors
from typing import Optional, Union

from pharmacophore.constants import is_hbond_acceptor_nitrogen, resolve_pml_output_paths


def _cliff_atom_roles(atom: Chem.Atom) -> tuple[bool, bool]:
    """Return (is_donor, is_acceptor) for a heavy cliff atom (handles explicit and implicit hydrogens)."""
    symbol = atom.GetSymbol()
    donor = symbol in ["N", "O", "S"] and atom.GetTotalNumHs(includeNeighbors=True) > 0
    acceptor = (symbol == "O" and atom.GetFormalCharge() <= 0) or (
        symbol == "N" and is_hbond_acceptor_nitrogen(atom)
    )
    return donor, acceptor


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
        self.mols: Optional[list[Chem.Mol]] = None
        self._mcs_mol: Optional[Chem.Mol] = None
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
        :param bits: int
            Fingerprint size in bits. This param is ignored if 'maccs' is used.
        :param radius: int
            Param only used for 'morgan' and 'featmorgan' fingerprints.
        :param type: str
            Available fingerprint type: 'morgan' (default), 'featmorgan', 'rdkit', 'atompair', 'torsion', 'maccs',
            'pattern', 'layered' or 'avalon'.
        :return: pd.DataFrame
        """
        if smi_col is None:
            smi_col = self.smi_col

        data = self.data
        smi_list = data[smi_col].tolist()
        pIC50_list = data['pIC50'].tolist()

        # add fp col
        mol_list = [Chem.MolFromSmiles(x) for x in smi_list]
        data['fp'] = [_calculate_fingerprint(x, fp_type=type, radius=radius, fpSize=bits) for x in mol_list]
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
        self._mcs_mol = mcs_mol

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

        if savepath:
            if SVG:
                with open(savepath, "w") as f:
                    f.write(img.data if hasattr(img, "data") else str(img))
            else:
                if hasattr(img, "save"):
                    img.save(savepath)
                else:
                    with open(savepath, "wb") as f:
                        f.write(img.data)

        return img

    def _diff_indices_for_mol(self, mol: Chem.Mol, mol_idx: Optional[int] = None) -> list[int]:
        """
        Return heavy-atom indices in ``mol`` that fall outside the shared MCS.

        When ``highlight_cliffs()`` has stored ``self._mcs_mol``, matches the MCS directly onto ``mol``
        so differences reflect ``mol``'s own atom ordering (including 3D ``.sdf`` poses with explicit Hs).
        Falls back to per-molecule ``self.atom_difference`` indexing when ``self._mcs_mol`` is not set.
        """
        if self._mcs_mol is not None:
            match_atoms = mol.GetSubstructMatch(self._mcs_mol)
            if not match_atoms and any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
                tagged = Chem.Mol(mol)
                for atom in tagged.GetAtoms():
                    atom.SetIntProp("_sar_full_idx", atom.GetIdx())
                mol_no_h = Chem.RemoveHs(tagged)
                no_h_match = mol_no_h.GetSubstructMatch(self._mcs_mol)
                if no_h_match:
                    match_atoms = tuple(
                        mol_no_h.GetAtomWithIdx(i).GetIntProp("_sar_full_idx")
                        for i in no_h_match
                    )
            if match_atoms:
                match_set = set(match_atoms)
                return [
                    atom.GetIdx()
                    for atom in mol.GetAtoms()
                    if atom.GetAtomicNum() > 1 and atom.GetIdx() not in match_set
                ]

        if self.atom_difference is None:
            raise ValueError(
                "No activity-cliff atom differences found. Call highlight_cliffs() before "
                "view_cliffs() or output_cliffs()."
            )
        if not self.atom_difference:
            return []

        if mol_idx is not None:
            if 0 <= mol_idx < len(self.atom_difference):
                raw = self.atom_difference[mol_idx]
            elif len(self.atom_difference) == 1:
                raw = self.atom_difference[0]
            else:
                raise IndexError(
                    f"mol_idx={mol_idx} out of range for atom_difference of length {len(self.atom_difference)}."
                )
        elif len(self.atom_difference) == 1:
            raw = self.atom_difference[0]
        else:
            raise IndexError(
                "Cannot infer molecule index for atom_difference when multiple entries exist."
            )

        num_atoms = mol.GetNumAtoms()
        return [idx for idx in raw if 0 <= idx < num_atoms]

    def _resolve_mol_names(
        self,
        mol_list: list[Chem.Mol],
        mol_names: Optional[Union[str, list[str], tuple[str, ...]]] = None,
    ) -> list[str]:
        """Resolve molecule names from explicit input, ``_Name``, ``self.data['name']``, or canonical SMILES."""
        if mol_names is not None:
            if isinstance(mol_names, str):
                if len(mol_list) != 1:
                    raise ValueError(
                        f"Expected {len(mol_list)} molecule names in mol_names, got single string {mol_names!r}."
                    )
                return [mol_names]
            if isinstance(mol_names, (list, tuple)):
                if len(mol_names) != len(mol_list):
                    raise ValueError(
                        f"Length of mol_names ({len(mol_names)}) must match number of molecules ({len(mol_list)})."
                    )
                return ["" if n is None else str(n) for n in mol_names]
            raise TypeError("mol_names must be a string or a sequence of strings.")

        has_name_col = "name" in self.data.columns
        smi_to_name: dict[str, str] = {}
        if has_name_col and self.smi_col in self.data.columns:
            for smi_val, name_val in zip(self.data[self.smi_col].tolist(), self.data["name"].tolist()):
                ref_mol = Chem.MolFromSmiles(str(smi_val)) if smi_val is not None else None
                if ref_mol is not None and pd.notna(name_val):
                    can_smi = Chem.MolToSmiles(ref_mol, canonical=True)
                    name_str = str(name_val).strip()
                    if name_str and can_smi not in smi_to_name:
                        smi_to_name[can_smi] = name_str

        resolved: list[str] = []
        for i, mol in enumerate(mol_list):
            prop_name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
            if prop_name:
                resolved.append(prop_name)
                continue
            if has_name_col and len(mol_list) == len(self.data):
                row_name = self.data["name"].iloc[i]
                if pd.notna(row_name) and str(row_name).strip():
                    resolved.append(str(row_name).strip())
                    continue
            mol_no_h = Chem.RemoveHs(Chem.Mol(mol))
            can_mol = Chem.MolToSmiles(mol_no_h, canonical=True)
            if can_mol in smi_to_name:
                resolved.append(smi_to_name[can_mol])
                continue
            resolved.append(f"Mol_{i}")
        return resolved

    def output_cliffs(
        self,
        mols: Optional[Union[Chem.Mol, list[Chem.Mol], tuple[Chem.Mol, ...]]] = None,
        savepath: Optional[Union[str, os.PathLike]] = None,
        mol_names: Optional[Union[str, list[str]]] = None,
    ) -> Union[str, list[str]]:
        """
        Output the 3D activity cliffs for one or more molecules as PyMOL ``.pml`` script(s).

        For a single molecule, ``savepath`` may be a ``.pml`` file path or an output directory
        (writing ``<dir>/<mol_name>_cliffs.pml``). For multiple molecules, ``savepath`` must be an
        output directory and writes ``<mol_name>_cliffs.pml`` per molecule with molecule-prefixed
        PyMOL object and group names.
        :param mols: Optional[Union[Chem.Mol, list[Chem.Mol], tuple[Chem.Mol, ...]]]
            Single 3D RDKit ``Chem.Mol``, sequence of ``Chem.Mol`` objects, or ``None`` to use the
            molecules previously passed to ``view_cliffs()``.
        :param savepath: Optional[Union[str, os.PathLike]]
            Required ``.pml`` file path (single molecule) or output directory path (single or multiple molecules).
        :param mol_names: Optional[Union[str, list[str]]]
            Optional molecule name(s) used for ``<mol_name>_cliffs.pml`` filenames and PyMOL object prefixes.
        :return: Union[str, list[str]]
            Written ``.pml`` path for a single molecule, or list of written paths for multiple molecules.
        """
        if mols is None:
            if self.mols is None:
                raise ValueError("No molecules available. Pass mols or call view_cliffs() first.")
            mol_list = list(self.mols)
        elif isinstance(mols, Chem.Mol):
            mol_list = [mols]
        elif isinstance(mols, (list, tuple)):
            if not mols:
                raise ValueError("mols sequence cannot be empty.")
            for idx, m in enumerate(mols):
                if not isinstance(m, Chem.Mol):
                    raise TypeError(f"Expected RDKit Chem.Mol at mols[{idx}], got {type(m).__name__}.")
            mol_list = list(mols)
        else:
            raise TypeError("mols must be an RDKit Chem.Mol, a sequence of Chem.Mol objects, or None.")

        if self._mcs_mol is None and self.atom_difference is None:
            raise ValueError(
                "No activity-cliff atom differences found. Call highlight_cliffs() before output_cliffs()."
            )

        for idx, m in enumerate(mol_list):
            if m.GetNumConformers() == 0:
                raise ValueError(
                    f"Molecule at index {idx} has no 3D conformer. Embed or load 3D coordinates before output_cliffs()."
                )

        names = self._resolve_mol_names(mol_list, mol_names=mol_names)
        out_paths, safe_names = resolve_pml_output_paths(savepath, names, suffix="cliffs")

        is_multi = len(out_paths) > 1
        for idx, (mol, out_path, mol_prefix) in enumerate(zip(mol_list, out_paths, safe_names)):
            diff_indices = self._diff_indices_for_mol(mol, mol_idx=idx)
            _write_single_cliff_pml(
                mol=mol,
                diff_indices=diff_indices,
                mol_prefix=mol_prefix,
                savepath=out_path,
                verbose=not is_multi,
            )

        if is_multi:
            print(f"Number of molecules: {len(out_paths)}")
            print(f"Activity cliff scripts written to {str(savepath).strip()}.")

        return out_paths if is_multi else out_paths[0]

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

    def _cliff_shapes(self, mol, mol_idx: Optional[int] = None):
        """
        Support function for _render_cliffs and _marimo_cliffs. Computes the highlight spheres and the per-atom sphere
        styles for the atoms of one molecule that fall outside the MCS.
        """
        if mol.GetNumConformers() == 0:
            raise ValueError(
                "Molecule has no 3D conformer. Embed or load 3D coordinates before rendering activity cliffs."
            )

        diff_indices = self._diff_indices_for_mol(mol, mol_idx=mol_idx)

        # track aromatic atoms
        spheres = []
        aromatic_atom_indices = set()
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            # check if aromatic atom not in MCS
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                if any(idx in diff_indices for idx in ring):
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
        for idx in diff_indices:
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            if symbol == 'H':
                continue

            # tag donor/acceptor
            donor, acceptor = _cliff_atom_roles(atom)

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

        spheres, atoms = self._cliff_shapes(mol, mol_idx=index)
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
            spheres, atoms = self._cliff_shapes(mol, mol_idx=i)
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


def _write_single_cliff_pml(
    mol: Chem.Mol,
    diff_indices: list[int],
    mol_prefix: str,
    savepath: str,
    verbose: bool = True,
) -> None:
    """Write a PyMOL ``.pml`` script highlighting one molecule's non-MCS activity-cliff features."""
    with open(savepath, "w") as f:
        f.write("# Drag-and-drop PyMOL script generated by pharmacophore-toolkit (SAR)\n")
        f.write("set auto_zoom, off\n\n")

        f.write("set_color Hydrophob, [46, 204, 113]\n")
        f.write("set_color HDonor, [33, 150, 243]\n")
        f.write("set_color HAcceptor, [244, 67, 54]\n")
        f.write("set_color DualH, [255, 0, 255]\n")
        f.write("set_color Aromatic, [255, 235, 59]\n\n")

        conf = mol.GetConformer()
        aromatic_atom_indices = set()
        group_members: list[str] = []
        ring_info = mol.GetRingInfo()
        for x, ring in enumerate(ring_info.AtomRings()):
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                if any(idx in diff_indices for idx in ring):
                    coords = [conf.GetAtomPosition(idx) for idx in ring]
                    centroid = np.mean(coords, axis=0)

                    obj_name = f"{mol_prefix}_Aromatic_{x}"
                    f.write(
                        f"pseudoatom {obj_name}, pos=[{centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f}]\n"
                    )
                    f.write(f"show spheres, {obj_name}\n")
                    f.write(f"hide nonbonded, {obj_name}\n")
                    f.write(f"set sphere_scale, 0.5, {obj_name}\n")
                    f.write(f"color Aromatic, {obj_name}\n")
                    f.write(f"set sphere_transparency, 0.0, {obj_name}\n\n")
                    group_members.append(obj_name)

                    for idx in ring:
                        aromatic_atom_indices.add(idx)

        for idx in diff_indices:
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            if symbol == 'H':
                continue
            pos = conf.GetAtomPosition(idx)

            donor, acceptor = _cliff_atom_roles(atom)

            color = None
            label = ""
            if donor and acceptor:
                color, label = "DualH", "Dual"
            elif donor:
                color, label = "HDonor", "Donor"
            elif acceptor:
                color, label = "HAcceptor", "Acceptor"
            elif idx not in aromatic_atom_indices and symbol in ['C', 'F', 'Cl', 'Br', 'I']:
                color, label = "Hydrophob", "Hydrophobic"

            if color:
                obj_name = f"{mol_prefix}_Pin_{label}_{idx}"
                f.write(f"pseudoatom {obj_name}, pos=[{pos.x:.3f}, {pos.y:.3f}, {pos.z:.3f}]\n")
                f.write(f"show spheres, {obj_name}\n")
                f.write(f"color {color}, {obj_name}\n")
                f.write(f"set sphere_transparency, 0.3, {obj_name}\n")
                scale = 0.7 if idx in aromatic_atom_indices else 1.0
                f.write(f"set sphere_scale, {scale}, {obj_name}\n\n")
                group_members.append(obj_name)

        group_name = f"{mol_prefix}_cliffs"
        if group_members:
            f.write(f"group {group_name}, {' '.join(group_members)}\n")
        f.write("set auto_zoom, on\n")
        if group_members:
            f.write(f"zoom {group_name}\n")
        else:
            f.write("zoom all\n")

    if verbose:
        print(f"Activity cliff script written to {savepath}.")


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


def _calculate_fingerprint(mol, fp_type='morgan', radius=2, fpSize=1024):
    """
    Calculates a bit-vector fingerprint using one of RDKit's algorithms.
    :param mol: Chem.Mol
        RDKit Mol object.
    :param fp_type: str
        One of 'morgan' (ECFP), 'featmorgan' (FCFP), 'rdkit', 'atompair', 'torsion', 'maccs', 'pattern',
        'layered' or 'avalon'.
    :param radius: int
        Only used for Morgan and FeatMorgan types.
    :param fpSize: int
        Size of the bit vector. Parma is ignored if 'maccs' is used.
    :return: ExplicitBitVect
    """
    fingerprints = {
        'morgan': lambda: rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fpSize).GetFingerprint(mol),
        'featmorgan': lambda: rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=fpSize,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()).GetFingerprint(mol),
        'rdkit': lambda: rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=fpSize).GetFingerprint(mol),
        'atompair': lambda: rdFingerprintGenerator.GetAtomPairGenerator(fpSize=fpSize).GetFingerprint(mol),
        'torsion': lambda: rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=fpSize).GetFingerprint(mol),
        'maccs': lambda: MACCSkeys.GenMACCSKeys(mol),
        'pattern': lambda: Chem.PatternFingerprint(mol, fpSize=fpSize),
        'layered': lambda: Chem.LayeredFingerprint(mol, fpSize=fpSize),
        'avalon': lambda: pyAvalonTools.GetAvalonFP(mol, nBits=fpSize),
    }

    if fp_type.lower() not in fingerprints:
        raise ValueError(f"Unknown fingerprint type: {fp_type}. Choose from {list(fingerprints)}")
    return fingerprints[fp_type.lower()]()


def _color_to_rgb(color_input):
    """Support function to convert Matplotlib colors to RGB for RDKit highlighting"""
    try:
        return mcolors.to_rgb(color_input)
    except ValueError:
        return f"Error: '{color_input}' is not a recognized color name or hex code."


if __name__ == "__main__":
    import doctest

    doctest.testmod()
