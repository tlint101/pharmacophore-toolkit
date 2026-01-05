# A Simple Pharmacophore-Toolkit

[![Pharmacophore-Toolkit Versions](https://img.shields.io/pypi/v/pharmacophore-toolkit.svg?label=Pharmacophore-Toolkit&color=blue)](https://pypi.org/project/pharmacophore-toolkit/)
[![Python Versions](https://img.shields.io/badge/python-3.9+-blue.svg?logo=python&logoColor=white)](https://pypi.org/project/pharmacophore-toolkit/)
[![Pharmacophore-Toolsets Test](https://github.com/tlint101/pharmacophore-toolkit/actions/workflows/tests.yml/badge.svg)](https://github.com/tlint101/pharmacophore-toolkit/actions/workflows/tests.yml)
[![Read the Docs](https://readthedocs.org/projects/pharmacophore-toolkit/badge/?version=latest)](https://pharmacophore-toolkit.readthedocs.io/en/latest/#)


The Pharmacophore-Toolkit is built on RDKit and allows for building simple pharmacophore models. The
Pharmacophore-Toolset can generate models from crystal structures, docking poses, or SMILES string. To generate a 3D
model, a .pml file will be generated. This files contains scripts to generate spheres with color and XYZ coordinates
defined. The final 3D image can be rendered in PyMOL. 

Documentation for the project can be found [here](https://pharmacophore-toolkit.readthedocs.io/en/latest/#).

## Install

You can install the Pharmacophore-Toolkit using pip:

```
pip install pharmacophore-toolkit
```

Alternatively, the environment can be created by cloning the repository and then running the following conda script:

```
conda env create -f environment.yaml
```

> [!NOTE]
> **Note:** The Pharmacophore-Toolkit relies on cairosvg to create images before being converted into .png format.
> Depending on your
> workstation/machine, the [CairoSVG](https://github.com/Kozea/CairoSVG) package will need to be installed manually.
> Installation instructions can be found [here](https://cairosvg.org). If it is not already installed on your machine
> globally,
> users can try conda to install cairosvg:

```
conda install conda-forge::cairosvg
```

Conda installation solved the problem on my machine, but may differ with yours.

## Tutorials

Tutorials are written as JupyterNotebooks and can be found [here](https://pharmacophore-toolkit.readthedocs.io/en/latest/tutorials/tutorials.html). The Pharmacophore-Toolkit can generate 
several types of images:

## Example Images

### 3D Model

A 3D conformation of molecules and alignment was performed using RDKit. Spheres can be generated for each query
molecule to highlight pharmacophore features, such as Hydrogen Bond Donor/Acceptor, Hydrophobic, Aromatic, etc. The
molecules can be generated in 3D using two methods:

- **[py3Dmol](https://github.com/3dmol/3Dmol.js)**

<figure>
    <img src="img/pharmacophore-demo.gif" width="400">
    <figcaption>
    The 3D model can be interactive in a Jupyter Notebook.
    </figcaption>
</figure>

<figure>
    <img src="img/3d_example_py3dmol.png" width="400">
    <figcaption>
    Images rendered in Jupyter Notebook using py3Dmol using a screenshot. Pharmacophore features for each molecule is 
    highlighted. Blue spheres represent Hydrogen Bond Donors, gold spheres for Aromatic rings, and green for 
    Hydrophobes. 
    </figcaption>
</figure>



- **[PyMOL](https://pymol.org)**

<figure>
    <img src="img/3d_example_pymol.png" width="400">
    <figcaption>
    Images rendered in PyMOL using the generated .pml file. Molecules are colored using default settings and mirror 
    those found using py3Dmol above. In this example, the spheres are mapped to Serotonin. Blue spheres represent 
    Hydrogen Bond Donors and gold spheres represent Aromatic rings.  
    </figcaption>
</figure>

### 2D Model

The pharmacophore features can also be depicted as a 2D image, where each atom is highlighted with their respective
pharmacophore features.

<figure>
    <img src="img/2d_pharmacophore.png" width="400">
    <figcaption>
    2D Pharmacophore images. Images were lightly edited to remove duplicate figure legends. 
    </figcaption>
</figure>

### Similarity Map

For funsies, a method to generate [Similarity Maps](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43)
are also available. Final image generation requires Cairo installation. Additional information on how Similarity Maps
work can be found on
the [RDKit blog](https://greglandrum.github.io/rdkit-blog/posts/2020-01-03-similarity-maps-with-new-drawing-code.html).

<figure>
    <img src="img/similarity_exmaple.png" width="400">
    <figcaption>
    Similarity map of molecules. All molecules were compared to Serotonin. More information can be seen
    at the <a href="https://www.rdkit.org/docs/source/rdkit.Chem.Draw.SimilarityMaps.html">RDKit documentation</a>.
    </figcaption>
</figure>