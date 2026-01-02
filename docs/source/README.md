# A Simple Pharmacophore-Toolkit

[![Pharmacophore-Toolkit Versions](https://img.shields.io/pypi/v/pharmacophore-toolkit.svg?label=Pharmacophore-Toolkit&color=blue)](https://pypi.org/project/pharmacophore-toolkit/)
[![Python Versions](https://img.shields.io/pypi/pyversions/pharmacophore-toolkit?style=flat&logo=python&logoColor=white)](https://pypi.org/project/pharmacophore-toolkit/)


The Pharmacophore-Toolkit is built on RDKit and allows for building simple pharmacophore models. The
Pharmacophore-Toolset can generate models from crystal structures, docking poses, or SMILES string. To generate a 3D
model, a .pml file will be generated. This files contains scripts to generate spheres with color and XYZ coordinates
defined. The final 3D image can be rendered in PyMOL. Additional information can be found under the
[Tutorials](tutorials/tutorials.md) section.

## Example Images

### 3D Model

A 3D conformation of molecules and alignment was performed using RDKit. Spheres can be generated for each query
molecule to highlight pharmacophore features, such as Hydrogen Bond Donor/Acceptor, Hydrophobic, Aromatic, etc. The
molecules can be generated in 3D using two methods:

**[py3Dmol](https://github.com/3dmol/3Dmol.js)**

```{figure} img/pharmacophore-demo.gif
:alt: py3Dmol
:width: 800px
:align: center
The 3D model can be interactive in a Jupyter Notebook.
```

```{figure} img/3d_example_py3dmol.png
:alt: py3Dmol
:width: 800px
:align: center
Images rendered in Jupyter Notebook using py3Dmol using a screenshot. Pharmacophore features for each molecule is 
    highlighted. Blue spheres represent Hydrogen Bond Donors, gold spheres for Aromatic rings, and green for 
    Hydrophobes. 
```

**[PyMOL](https://pymol.org)**

```{figure} img/3d_example_pymol.png
:alt: PyMOL
:width: 800px
:align: center
Images rendered in PyMOL using the generated .pml file. Molecules are colored using default settings and mirror 
    those found using py3Dmol above. In this example, the spheres are mapped to Serotonin. Blue spheres represent 
    Hydrogen Bond Donors and gold spheres represent Aromatic rings.  
```

### 2D Model

The pharmacophore features can also be depicted as a 2D image, where each atom is highlighted with their respective
pharmacophore features.

```{figure} img/2d_pharmacophore.png
:alt: 2d_pharmacophore
:width: 800px
:align: center
2D Pharmacophore images. Images were lightly edited to remove duplicate figure legends. 
```

### Similarity Map

For funsies, a method to generate [Similarity Maps](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43)
are also available. Final image generation requires Cairo installation. Additional information on how Similarity Maps
work can be found on
the [RDKit blog](https://greglandrum.github.io/rdkit-blog/posts/2020-01-03-similarity-maps-with-new-drawing-code.html).

```{figure} img/similarity_exmaple.png
:alt: similarity_map
:width: 400px
:align: center
Similarity map of molecules. All molecules were compared to Serotonin. More information can be seen
    at the <a href="https://www.rdkit.org/docs/source/rdkit.Chem.Draw.SimilarityMaps.html">RDKit documentation</a>.
```
