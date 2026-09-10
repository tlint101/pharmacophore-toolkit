import marimo

__generated_with = "0.24.0"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Marimo - Activity Cliffs
    The following notebook demos usage of the pharmacophore-toolkit activity cliffs in Marimo notebooks. This notebook essentially follows the same format and scripts as seen in the Jupyter notebook tutorial [03-activity_cliffs.ipynb]([03-.ipynb](https://github.com/tlint101/pharmacophore-toolkit/blob/main/tutorials/03-activity_cliffs.ipynb))

    **NOTE:** This notebook will not render in GitHub. It should be downloaded on your machine and rendered in Marimo instead.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem.Draw import MolsToGridImage
    from rdkit.Chem import AllChem, rdMolAlign
    from pharmacophore import SAR

    return AllChem, Chem, MolsToGridImage, SAR, mo, pd, rdMolAlign


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Prepare Dataset

    Molecules for this demo will be [Zaprinast](https://en.wikipedia.org/wiki/Zaprinast) and [Sildenafil](https://en.wikipedia.org/wiki/Sildenafil). The latter is more commonly known on teh market as Viagra, while the former is an unsuccessful clinical drug candidate and precursor. Both molecules are Phosphodiesterase 5 (PDE5) inhibitors. Their IC50 values were obtained from [MedChemExpress](https://www.medchemexpress.com).

    As the two are structurally related and contain a wide difference in pIC50 values (~2 log difference), they serve as an interesting test for the Structure Activity Relationship (SAR) module.
    """)
    return


@app.cell
def _(pd):
    data = {
        "name": ["zaprinastat", "sildenafil"],
        "smiles": [
            "O=C1C2=C(N=NN2)NC(C3=CC=CC=C3OCCC)=N1",
            "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C"
        ],
        "activity": ["1,220", "5.22"],
        "units": ["nM", "nM"]
    }

    df = pd.DataFrame(data=data)
    df
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Visualize Molecules
    Before we proceed, we will first visualize the molecules used in this demo. A quick visual inspection shows similar core scaffold between the two molecules. These differences could account for their difference potencies in PED5 inhibition.
    """)
    return


@app.cell
def _(Chem, MolsToGridImage, df):
    smi_list = df.smiles.tolist()
    mols = [Chem.MolFromSmiles(x) for x in smi_list]
    name_list = df.name.tolist()
    activity_list = df.activity.tolist()
    unit_list = df.units.tolist()
    label = []
    for x, y, z in zip(name_list, activity_list, unit_list):
        label.append(f"{x}\n{y}{z}")

    MolsToGridImage(mols=mols, legends=label, subImgSize=(300,300))
    return label, mols


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Structure Activity Relationship

    To explore the different structures of the two molecules, the **SAR()** class will be instantiated. The **SAR()** class requires a pd.Dataframe of the data with a columns for the structure smiles and activity. The activity is, presumably, the IC50 value for the query moelcule. Internally, the **SAR()** class will convert the activity column into the IC50 range.
    """)
    return


@app.cell
def _(SAR, df):
    # activity needs to be a float
    sar = SAR(data=df, smi_col='smiles', act_col='activity', type="marimo")
    sar.data
    return (sar,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Calculating LiPE
    Once the **SAR()** class is instantiated, the LiPE for each molecule can be calculated. This is done using the **calc_LiPE()** method and is done internally, appending the LiPE score to the DataFrame.
    """)
    return


@app.cell
def _(sar):
    lipe = sar.calc_LiPE()
    lipe
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Calculating SALI
    Additionally, the SALI score between a part of molecules can be Calculated. This will continue to use the same DataFrame as given and modified with the **calc_LiPE()** method. The **get_sali()** method will reformat the DataFrame. In this demo, the two molecules are combined into a single row containing the relevant information as well as teh calculated Tanimoto Similarity score and SALI score.
    """)
    return


@app.cell
def _(sar):
    sali = sar.get_sali()
    sali
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Visualizing Activity Cliffs - 2D
    The activity cliffs for a pair of molecules can be drawn using the **highlight_cliffs()** method. While a list of smiles can be given, if None is given, then the smiles column given at **SAR()** instantiation will be used.

    Internally, the **highlight_cliffs** will identify the Maximum Common Scaffolds (MCS) between the molecules. The different atoms will be identified and highlighted when the molecules are drawn.

    Alignment of the 2D molecules will be handled. However, for the demo molecules they do not share an exact core scaffolds. Thus, they are not properly aligned for this example.
    """)
    return


@app.cell
def _(sar):
    sar.highlight_cliffs(highlight_color="deepskyblue")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Visualizing Activity Cliffs - 3D

    The molecules and their activity cliffs can also be visualized in 3D. For this demo, conformations for the molecules will be generated.

    **Note:** The 3D conformations are only needed for this demo. However, it is not strictly needed for 3D visualization. Docking poses or other methods can be used as inputs as well.
    """)
    return


@app.cell
def _(AllChem, Chem, MolsToGridImage, label, mols, rdMolAlign):
    mols_h = [Chem.AddHs(m) for m in mols]
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    for m in mols_h:
        AllChem.EmbedMolecule(m,ps)

    aligned = []

    for a in mols_h:
        mol_aligned = rdMolAlign.GetO3A(a,mols_h[0])
        aligned.append(mol_aligned.Align())

     # remove hydrogens for clarity.
    mols_noH = [Chem.RemoveHs(m) for m in mols_h]
    MolsToGridImage(mols=mols_noH, legends=label, subImgSize=(300,300))
    return (mols_noH,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Visualize Activity Cliffs - py3DMol
    The molecules can be visualized in 3D using py3DMol. this is done using the **view_cliffs()** method. The molecules will have spheres highlighting the differences as seen with the **highlight_cliffs()** method. The colors correspond to:
    - HAcceptor -> Red
    - HDonor -> Blue
    - Both HAcceptor/Donor -> Magenta
    - Hydrophobic -> Green
    - Aromatic -> Gold
    """)
    return


@app.cell
def _(mols_noH, sar):
    sar.view_cliffs(mols_noH)
    return


if __name__ == "__main__":
    app.run()
