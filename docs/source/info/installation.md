# Installation
The recommended way to install Pharmacophore-Toolkit is through PyPI using pip:

```
pip install pharmacophore-toolkit
```

Alternatively, the environment can be created by cloning the repository and then running the following conda script:

```
conda env create -f environment.yaml
```

## Note

The Pharmacophore-Toolkit relies on 'cairosvg' to create images before being converted into '.png' format. Depending on
your workstation/machine, the '[CairoSVG](https://github.com/Kozea/CairoSVG)' may need to be installed.

Installation instructions can be found at [CairoSVG documentation](https://cairosvg.org).

If `cairosvg` is not installed globally, users can try installing it via `conda`:
```
conda install conda-forge::cairosvg
```
or via `pip`
```
pip install cairosvg
```
Conda installation solved the problem on my machine, but may differ with yours.
