from rdkit import Chem
from rdkit.Chem import AllChem
from pharmacophore import Pharmacophore

# script adapted from tutorial to check if it installs and runs

molecules = {"serotonin": "C1=CC2=C(C=C1O)C(=CN2)CCN",
             "psilocin": "CN(C)CCc1c[nH]c2cccc(O)c12",
             "mescaline": "O(c1cc(cc(OC)c1OC)CCN)C"}

def test():
    mol_smi = [x for x in molecules.values()]
    mol_name = [x for x in molecules.keys()]
    mols = [Chem.MolFromSmiles(x) for x in mol_smi]

    # generate conformations
    mols = [Chem.AddHs(m) for m in mols]
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    for m in mols:
        AllChem.EmbedMolecule(m,ps)

    mols_noH = [Chem.RemoveHs(m) for m in mols]
    pharm = Pharmacophore()
    df = pharm.to_df(mols_noH, mol_name)

    print(df)

if __name__ == "__main__":
    test()