"""
Script to draw pharmacophore
"""

import os
import io
import py3Dmol
import json
import warnings
import matplotlib.image as img
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from typing import Optional, Union, Any
from PIL import Image
from collections import defaultdict
from cairosvg import svg2png
from IPython.display import SVG
from pharmacophore.constants import FEATURE_COLORS, INTERACTIVE_COLORS, color_convert
from pharmacophore import Pharmacophore
from rdkit import Chem
from rdkit.Chem import rdDepictor, AllChem, Mol
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.Chem.Draw.MolDrawing import DrawingOptions


class Draw:
    def __init__(self, mol: Optional[Chem.Mol] = None, features: Optional[Union[str, dict]] = 'default'):
        self.mol = mol
        self.features = features

    def draw_pharm(self, mol: Chem.Mol, features: Optional[Union[str, dict]] = None, color: dict = None,
                   savepath: str = None):
        """
        Draw a 2D pharmacophore image
        :param mol: Chem.Mol
            A molecule in ROMol format.
        :param features: str
            Determines which feature algorithm to use. Can use 'default' or 'rdkit'.
        :param color: dict
            Options to modify color of features. Must contain the name of the features. Colors can be given as either
            distinct name or hex codes.
        :param savepath: str
            Set the savepath to save the image.
        :return:
        """
        # set instance variables
        global color_palette
        if mol is None:
            mol = self.mol
        if features is None:
            features = self.features

        # calculate pharmacophore features
        try:
            pharm = Pharmacophore(features=features)
            calc_features = pharm.calc_pharm(mol=mol)
        except:
            raise ValueError('Unknown features! Only "default" or "rdkit", or custom features as dict are supported!')

        # dictionaries to highlight atoms and highlight radius
        atom_highlights = defaultdict(list)
        highlight_rads = {}

        # set to default or custom colors
        if color is None:
            color_palette = FEATURE_COLORS
        elif isinstance(color, dict):
            for key in color:
                color[key] = color_convert(color[key])
            color_palette = color

        # calc features and append results to atom_highlights and highlight_rads
        for feature in calc_features:
            if feature[0] in color_palette:
                # extract colors from constants
                color = color_palette[feature[0]]
                # # for troubleshooting
                # print(color)
                for atom_id in feature[1]:
                    atom_highlights[atom_id].append(color)
                    highlight_rads[atom_id] = 0.5

        # flatten molecule into 2D
        rdDepictor.Compute2DCoords(mol)
        rdDepictor.SetPreferCoordGen(True)
        drawer = rdMolDraw2D.MolDraw2DSVG(800, 800)

        # set drawing options
        # Use black for all elements
        drawer.drawOptions().updateAtomPalette(
            {k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()}
        )
        drawer.SetLineWidth(2)
        drawer.SetFontSize(6.0)
        drawer.drawOptions().continuousHighlight = False
        drawer.drawOptions().splitBonds = False
        drawer.drawOptions().fillHighlights = True

        # get atom label
        for atom in mol.GetAtoms():
            atom.SetProp("atomLabel", atom.GetSymbol())
        drawer.DrawMoleculeWithHighlights(
            mol, "", dict(atom_highlights), {}, highlight_rads, {}
        )
        drawer.FinishDrawing()

        # draw molecule and save a temporary file
        svg = drawer.GetDrawingText().replace("svg:", "")
        SVG(svg)
        with open(f"pharm.svg", "w") as f:
            f.write(svg)

        # convert svg into png
        svg2png(bytestring=svg, write_to=f"image.png")

        # Set figure legend for feature and color type
        fig, (ax, picture) = plt.subplots(
            nrows=2,
            figsize=(4, 4),
            gridspec_kw={"height_ratios": [1, 5]},
        )

        # draw image and remove temporary file
        mol_image = img.imread(f"image.png")
        picture.imshow(mol_image)
        picture.axis("off")
        os.remove(f"image.png")
        os.remove(f"pharm.svg")

        # Data for the circles
        circle_radii = [0, 50, 100, 150, 200, 250]
        feature_values = list(color_palette.values())  # extract color values
        circle_colors = [i for i in feature_values]  # match circle to color values

        # get color labels for fig legend
        if isinstance(features, dict):
            circle_annotations = list(features.keys())
        elif features == "rdkit":
            circle_annotations = [
                "Donor",
                "Acceptor",
                "Aromatic",
                "Hydrophobe",
                "LumpedHydrophobe",
                "PosIonizable",
            ]
        else:
            circle_annotations = [
                "Donor",
                "Acceptor",
                "Aromatic",
                "Hydrophobe",
            ]
        # Draw the circles and annotations
        fontsize = 6
        if features == "rdkit":
            fontsize = 4

        for radius, color, annotation in zip(
                circle_radii, circle_colors, circle_annotations
        ):
            x = radius
            circle = plt.Circle((x, -5), 5, color=color)  # , alpha=0.5)
            ax.add_patch(circle)
            ax.annotate(
                annotation,
                (x, 10),
                va="center",
                ha="center",
                fontsize=fontsize,
                fontweight="bold",
            )

        # Set axis limits
        ax.set_xlim(-10, 270)
        ax.set_ylim(-20, 20)
        ax.axis("off")

        # Set aspect ratio to equal
        ax.set_aspect("equal", adjustable="box")
        if savepath:
            plt.savefig(f"{savepath}", dpi=300)

    # support function to draw molecule with atom index
    def atom_number(self, mol: Optional[Chem.Mol] = None, label: str = "atomNote", size: tuple = (300, 300)):
        """
        Draw query molecule with labeled RDKit atom indices.
        :param mol: Chem.Mol
            A molecule in ROMol format.
        :param label: str
            Determines which style to label the molecule. Defaults to atomNote. In total, can use 'atomNote',
            'atomLabel', and 'molAtomMapNumber'.
        :param size: tuple
            Determine the size of the molecule to draw.
        :return:
        """
        # set instance variable
        if mol is None:
            mol = self.mol

        # Check if 2D coordinates are missing, and compute them if necessary. This will also strip away 3D coordinates.
        if not mol.GetNumConformers() or mol.GetConformer().Is3D():
            AllChem.Compute2DCoords(mol)

        # get atom indices
        for atom in mol.GetAtoms():
            if label == 'atomNote' or label == 'atomLabel' or label == 'molAtomMapNumber':
                atom.SetProp(label, str(atom.GetIdx()))
            else:
                raise ValueError("Only 'atomNote', 'atomLabel', and 'molAtomMapNumber' accepted!")

        # draw molecule
        img = Chem.Draw.MolToImage(mol, size=size)

        return img

    def similarity_maps(self, refmol: Chem.Mol = None, querymol: Chem.Mol = None, radius: int = 2, nbits: int = 2048,
                        fpType: str = 'bv', cmap: Optional[Union[str, list]] = None, savepath: str = None):
        """
        Generate similarity map between query molecule and reference molecule. Only Morgan fingerprints are used.
        Image must be in png format.
        :param refmol: Chem.Mol
            Reference molecule. Must be in ROMol format.
        :param querymol: Chem.Mol
            Query molecule. Must be in ROMol format.
        :param radius: int
            Set the radius for the MorganFingerprint. Default to 2.
        :param nbits: int
            Set the number of bits for the generated fingerprint. Default to 2048.
        :param fpType: str
            Set the fingerprint type. Default to 'bv'. Can only use 'bv' or 'count'
        :param cmap: Optional[Union[str, list]]
            Set the colormap of figure. Uses matplotlib color scheme: https://matplotlib.org/stable/users/explain/colors/colormaps.html
        :param savepath: str
            Set image savepath.
        :return:
        """
        # set instance variable
        if refmol is None:
            refmol = self.mol

        # check fpType
        if fpType not in ['bv', 'count']:
            raise ValueError("Only 'bv' or 'count' accepted!")

        # convert list into a matplotlib colormap
        if isinstance(cmap, list):
            if len(cmap) != 3:
                raise ValueError("Only three colors are accepted!")
            cmap = LinearSegmentedColormap.from_list("custom_similarity_colors", cmap)

        d = Chem.Draw.MolDraw2DCairo(400, 400)
        function = lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=radius, fpType=fpType, nBits=nbits)
        _, maxWeight = SimilarityMaps.GetSimilarityMapForFingerprint(refMol=refmol, probeMol=querymol,
                                                                     fpFunction=function, draw2d=d, colorMap=cmap)

        # finish drawing
        d.FinishDrawing()
        img_data = d.GetDrawingText()

        # draw image as png file
        bio = io.BytesIO(img_data)
        img = Image.open(bio)

        # if savepath given in params
        if savepath:
            img.save(savepath)

        return img


class View:
    def __init__(self, mol: Optional[Union[Chem.Mol, list[Chem.Mol]]] = None,
                 pharmacophore: Optional[Union[str, dict]] = 'default', type: str = "jupyter"):
        """
        Shared functions for the View() class.
        :param mol: Optional[Union[Chem.Mol, list[Chem.Mol]]]
            A set of ROMols in a list to render.
        :param pharmacophore: Optional[Union[str, dict]]
            A set of pharmacophores for rendering.
        :param type: str
            Set the output for Jupyter or Marimo notebooks. Defaults to Jupyter.
        """
        # check type
        if type.lower() not in ("jupyter", "marimo"):
            raise ValueError("Only 'jupyter' or 'marimo' accepted!")
        else:
            self.type = type
            self.mol = mol
            self.pharmacophore = pharmacophore

    def view(self, mol: Optional[Union[list[Chem.Mol], Chem.Mol]] = None, pharmacophore: Optional[list] = None,
             color: Optional[dict] = None, labels: bool = True, window: tuple = (500, 500),
             prefix: str = "pharmacophore"):
        """
        Generate an interactive py3Dmol image in either Jupyter or Marimo Notebook of the molecule and its
        pharmacophores. Must include a list of molecules and a list of pharmacophores generated using
        Pharmacophore.calc_pharm().
        :param mol: Optional[Union[list[Chem.Mol], Chem.Mol]]
            A molecule in ROMol format. Can be a single molecule or a list of molecules.
        :param pharmacophore: Optional[list]
            A list containing a list of pharmacophore data generated from Pharmacophore.calc_pharm().
        :param color: Optional[dict]
            A dictionary containing the following: {pharmacophore:color}. The color name can be in hex code or color
            name. If None given, will use default colors.
        :param labels: bool
            Whether to generate labels overlay on the pharmacophore spheres.
        :param window: tuple
            Set the window size of the py3dmol figure.
        :param prefix: str
            Set the prefix for the saved image. Only works for Marimo notebooks.
        :return:
        """
        if self.type == "jupyter":
            return self._jupyter(mol, labels, color, pharmacophore, window)
        elif self.type == "marimo":
            return self._marimo(mol, labels, color, pharmacophore, window, prefix)
        else:
            raise ValueError("Only 'jupyter' or 'marimo' accepted!")

    def _jupyter(self, mol: Optional[list[Mol]], labels: bool, color: Optional[dict[Any, Any]],
                 pharmacophore: Optional[list[Any]], window: tuple):
        # instantiate variables
        if pharmacophore is None:
            pharmacophore = self.pharmacophore

        # set to default or custom colors
        if color is None:
            color = INTERACTIVE_COLORS
        elif isinstance(color, dict):
            for key in color:
                color[key] = color_convert(color[key])
            color = color

        # convert rdkit mol into mol_block
        if not isinstance(mol, list):
            mol = [mol]

        # set vars to self for use in _render
        self.mols = mol
        self.pharmacophore = pharmacophore
        self.colors = color
        self.labels = labels
        self.window = window

        # dropdown menu
        import ipywidgets as widgets
        dropdown = widgets.Dropdown(
            options=[(f"Molecule {i + 1}", i) for i in range(len(mol))],
            value=0,
            description="Select:",
            style={"description_width": "initial"}
        )

        widgets.interact(self._render, index=dropdown)

    def _render(self, index):
        """Function to render molecules and its pharmacophore"""
        mol = self.mols[index]
        mol_block = Chem.MolToMolBlock(mol)

        viewer = py3Dmol.view(width=self.window[0], height=self.window[1])
        viewer.setBackgroundColor("white")
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()

        # map pharmacophore based on list of list inputs
        if isinstance(self.pharmacophore[0][0], str):
            pharma_index = self.pharmacophore  # if multiple pharmacophore is given
        else:
            pharma_index = self.pharmacophore[index]  # single pharmacophore
        for pharma in pharma_index:
            label = pharma[0]
            x = pharma[2]
            y = pharma[3]
            z = pharma[4]
            color = self.colors.get(label, )

            viewer.addSphere({
                "center": {"x": x, "y": y, "z": z},
                "radius": 0.5,
                "color": color,
                "opacity": 0.9,
            })

            if self.labels:
                viewer.addLabel(label, {
                    "position": {"x": x, "y": y, "z": z},
                    "fontSize": 12,
                    "showBackground": True,
                    "backgroundColor": color,
                    "fontColor": "white"
                })

        viewer.show()

    def _marimo(self, mol: Optional[list[Mol]], labels: bool, color: Optional[dict[Any, Any]],
                pharmacophore: Optional[list[Any]], window: tuple, prefix):
        """
        Output Marimo interactive window.
        """
        try:
            import marimo as mo
        except Exception as e:
            return e
        if not isinstance(mol, list):
            mol = [mol]
        if pharmacophore is None:
            calc = Pharmacophore(features=self.pharmacophore)
            pharmacophore = [calc.calc_pharm(mol=m) for m in mol]
        if color is None:
            color = INTERACTIVE_COLORS
        elif isinstance(color, dict):
            color = {key: color_convert(value) for key, value in color.items()}

        # flatten feature list
        if pharmacophore and pharmacophore[0] and isinstance(pharmacophore[0][0], str):
            pharma = [pharmacophore] * len(mol)
        else:
            pharma = pharmacophore

        # serialize molecule data
        data = []
        for i, (m, feats) in enumerate(zip(mol, pharma)):
            if m.GetNumConformers() == 0:
                warnings.warn(
                    f"Molecule {i + 1} has no conformer. calc_pharm() cannot compute "
                    f"feature centroids without 3D coordinates, so no spheres will be "
                    f"drawn. Embed the molecule first (AllChem.EmbedMolecule).",
                    stacklevel=2,
                )
            elif not feats:
                warnings.warn(f"Molecule {i + 1} has no pharmacophore features.", stacklevel=2)

            data.append({
                "name": f"Molecule {i + 1}",
                "mol": Chem.MolToMolBlock(m),
                "feats": [{"label": f[0],
                           "x": float(f[2]), "y": float(f[3]), "z": float(f[4]),
                           "color": self._css_color(color.get(f[0]))}
                          for f in feats],
            })

        width, height = window

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
        const LABELS = {json.dumps(bool(labels))};
        const box = document.getElementById("viewer");
        const sel = document.getElementById("sel");
        DATA.forEach((d, i) => sel.add(new Option(d.name, i)));
        const viewer = $3Dmol.createViewer(box, {{backgroundColor: "white"}});

        function render(i) {{
          const d = DATA[i];
          viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllLabels();
          viewer.addModel(d.mol, "mol");
          viewer.setStyle({{}}, {{stick: {{}}}});
          d.feats.forEach(f => {{
            const pos = {{x: f.x, y: f.y, z: f.z}};
            viewer.addSphere({{center: pos, radius: 0.5, color: f.color, opacity: 0.9}});
            if (LABELS) {{
              viewer.addLabel(f.label, {{position: pos, fontSize: 12, showBackground: true,
                                         backgroundColor: f.color, fontColor: "white"}});
            }}
          }});
          viewer.zoomTo(); viewer.render();
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

    @staticmethod
    def _css_color(c):
        """
        Convert color for css.
        Accept 'royalblue', '#4169e1' or an (r, g, b) float tuple from color_convert.
        """
        if isinstance(c, (tuple, list)):
            return "#%02x%02x%02x" % tuple(int(round(255 * v)) for v in c[:3])
        return c


if __name__ == "__main__":
    import doctest

    doctest.testmod()
