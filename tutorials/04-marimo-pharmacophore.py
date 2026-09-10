import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Pharmacophore-Toolkit in Marimo
    The following notebook demos usage of the pharmacophore-toolkit in Marimo notebooks. This notebook essentially follows the same format and scripts as seen in the Jupyter notebook tutorial [01-pharmacophore_tutorial.ipynb]([01-pharmacophore_tutorial.ipynb](https://github.com/tlint101/pharmacophore-toolkit/blob/main/tutorials/01-pharmacophore_tutorial.ipynb))

    **NOTE:** This notebook will not render in GitHub. It should be downloaded on your machine and rendered in Marimo instead.
    """)
    return


@app.cell
def _():
    import marimo as mo
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    from rdkit.Chem.Draw import MolsToGridImage
    from pharmacophore import Pharmacophore, Draw, View
    from pharmacophore.constants import INTERACTIVE_COLORS
    import py3Dmol

    return (
        AllChem,
        Chem,
        Draw,
        MolsToGridImage,
        Pharmacophore,
        View,
        mo,
        rdMolAlign,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Dataset

    This tutorial will use serotonin, psilocin, and mescaline as examples. All three binds to the Serotonin 5HT-2A receptor.

    In this example, the molecules will be given as a dictionary in a name:SMILES pair. Each pair is then extracted into a separate list containing the SMILES string and the molecule name. The molecules can be visualized for confirmation.
    """)
    return


@app.cell
def _(Chem, MolsToGridImage):
    molecules = {"serotonin": "C1=CC2=C(C=C1O)C(=CN2)CCN",
                 "psilocin": "CN(C)CCc1c[nH]c2cccc(O)c12",
                 "mescaline": "O(c1cc(cc(OC)c1OC)CCN)C"}

    mol_smi = [x for x in molecules.values()]
    mol_name = [x for x in molecules.keys()]
    mols = [Chem.MolFromSmiles(x) for x in mol_smi]

    MolsToGridImage(mols=mols, legends=mol_name)
    return mol_name, mols


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Generating Conformations
    The 3D conformation can now be generated. For this example, a randomseed is given to ensure reproducibility.
    """)
    return


@app.cell
def _(AllChem, Chem, mols):
    mols_h = [Chem.AddHs(m) for m in mols]
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    for m in mols_h:
        AllChem.EmbedMolecule(m,ps)
    return (mols_h,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Once the 3D conformations have been generated, the molecules can be aligned. This is done using RDKit's Open3D align implementation. The alignment results can also output the RMSD score between the alignments. In this case, the serotonin molecule is used as the "base" for alignment. The resulting RMSD score gives an indication of the alignment results.
    """)
    return


@app.cell
def _(mols_h, rdMolAlign):
    aligned = []

    ref = mols_h[0]

    for x in mols_h:
        mol_aligned = rdMolAlign.GetO3A(x, ref)
        aligned.append(mol_aligned.Align())

    # print(aligned)  # Optional. Will display a list of RMSD alignment score.
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### Optional
    Because the molecules in this tutorial were rendered from SMILES strings, the 3D conformation will need to be saved. In this case, the conformation will be saved in .sdf format. This is done as follows:
    """)
    return


@app.cell
def _(Chem, mol_name, mols_h):
    for mol, name in zip(mols_h, mol_name):
        w = Chem.SDWriter(f"data/{name}.sdf")
        w.write(mol)
        w.close()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Checking Conformation
    Molecule conformation can be checked from the generated .sdf script above. Or the molecule conformation can be generated in 2D.
    """)
    return


@app.cell
def _(Chem, MolsToGridImage, mol_name, mols_h):
    mols_noH = [Chem.RemoveHs(m) for m in mols_h]  # remove hydrogens for clarity.
    MolsToGridImage(mols=mols_noH, legends=mol_name)
    return (mols_noH,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Count Pharmacophores
    Pharmacophores for each molecule can be generated and displayed as a DataFrame. To generate the table, initialize the Pharmacophore() class. By default, the Pharmacophore-Toolkit features will be set to use "default" methods. The Pharmacophore-Toolkit can also utilize RDKit feature sets by setting the feature param to "rdkit" or a custom feature set given in a dictionary format.

    Once initialized, the to_df() method can be used. This method requires a list of molecules in Chem.Mol format and a list of the molecule name. The to_df() method also contains a "features" param that, if given, will override the param given upon Pharmacophore initialization. Again, "default", "rdkit", or a custom feature set in a dictionary format can be given. The generated table will display the number of matches to the given features in a DataFrame.
    """)
    return


@app.cell
def _(Pharmacophore, mol_name, mols_noH):
    pharma = Pharmacophore()

    df = pharma.to_df(mols_noH, mol_name)
    df
    return (pharma,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### Optional—Custom Features
    If the "default" or "RDKit" features are not suitable for a project, users can set their own features. This requires a dictionary in 'feature_name':[SMARTS string] format. An example is given to identify all aromatic carbon atoms.
    """)
    return


@app.cell
def _(mol_name, mols_noH, pharma):
    custom_feat = {'aromatic': ['c']}

    custom_df = pharma.to_df(mols_noH, mol_name, features=custom_feat)
    custom_df
    return (custom_feat,)


@app.cell
def _(pharma):
    pharma.features
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Checking Features
    Features used for a given model can be dispalyed using the feature_types() method.
    Output default features. Can be used to output features associated with rdkit using the "features" param.

    **NOTE:** Setting the feature types will overwrite the stored value. In this case, it is set to the custom_feat dictionary from above. In this example, the Pharmacophore() class will be reinitialized to use the default features.
    """)
    return


@app.cell
def _(Pharmacophore):
    pharm = Pharmacophore()
    print(pharm.feature_types())
    return (pharm,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Generate Pharmacophore Features
    The following script will generate pharmacophore features based on a query molecule. This will take the molecule with the 3D conformaiton generated above and output a list of lists, containing the feature type, matching atom index, and the XYZ position. The features for the serotonin molecule is given as an example.

    Again, the calc_pharm contains a "features" param that can modify the type of features for a molecule. As none is given at for this method or at the initialization of the Pharmacophore class, it will utilize default features.
    """)
    return


@app.cell
def _(mols_noH, pharm):
    # default
    pharmac = pharm.calc_pharm(mols_noH[0])
    pharmac
    return (pharmac,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Render Pharmacophore Features

    The pharmacophore-toolkit contains a View() class to quickly render the molecules and their pharmacophores in a Jupyter notebook. This requires calling the View() class and the view() method. The parameters includes an RDKit molecule or a list of RDKit molecules and the calculated pharmacophores genereated from calc_pharm().

    Additional parameters includes lables and window to label the pharmacophores and adjust the py3dmol window size, respectively.

    **NOTE:** This method defaults for Jupyter Notebook. To set to marimo notebooks, use arguemnt type.
    """)
    return


@app.cell
def _(View, mols_noH, pharmac):
    v = View(type="marimo")
    single = mols_noH[0]
    v.view(mols_noH, pharmac, labels=True, window=(500, 500))
    return (v,)


@app.cell
def _(mols_noH, pharm):
    pharma_list = []
    for y in mols_noH:
        calc = pharm.calc_pharm(y)
        pharma_list.append(calc)
    # pharma # To check the list of pharmacophores
    return (pharma_list,)


@app.cell
def _(mols_noH, pharma_list, v):
    v.view(mols_noH, pharma_list, labels=True, window=(500, 500))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Output Pharmacophore Features

    As of now, py3Dmol does not include a screenshot feature. Rendering the pharmacophores for each molecule may be easier in PyMOL. The calculated pharmacophore features can be saved as a .pml file and then opened in PyMOL directly with the molecules (in .sdf or .mol2 formate). This will allow users to render the images in higher quality and save the 3D orientation based on user preferences. In this example, the pharmacophores for mescaline will be generated.
    """)
    return


@app.cell
def _(pharm):
    pharm.output_features(savepath='data/pharma.pml')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Generate 2D Image
    Pharmacophores can also be highlighted on the molecule in 2D. This is done using the Draw() class. Much of the parameters are the same as those used in the Pharmacophore() class. The query molecule will be drawn and atoms matching the pharmacophore sets will be highlighted in the corresponding color.
    """)
    return


@app.cell
def _(Draw, mols_noH):
    drawing = Draw()
    drawing.draw_pharm(mols_noH[0])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    #### Optional Color Options
    Additionally, the colors on the 2D images can be modified by providing a dictionary for the color parameter in draw_pharm(). The dictionary should be in feature:color format. The feature must match the style given. i.e. same capitalization, spelling, etc. In this example, the custom_feat contains "aromatic" in lower case.
    """)
    return


@app.cell
def _(Draw, custom_feat, mols_noH):
    custom_draw = Draw(features=custom_feat)
    colors = {"aromatic":'gold'}
    custom_draw.draw_pharm(mols_noH[0], color=colors)
    return


if __name__ == "__main__":
    app.run()
