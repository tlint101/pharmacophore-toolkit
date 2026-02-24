import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Crippen import MolLogP
from typing import Optional


class SAR:
    def __init__(self, data: pd.DataFrame, smi_col: str = "smiles", act_col: str = "activity", units: str = "nM"):
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
        LipE = pIC50 - LogP
        https://en.wikipedia.org/wiki/Lipophilic_efficiency
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

    def get_sali(self, smi_col: Optional[list] = None, act_col: Optional[list] = None, ic50: bool = False,
                 units: str = "nM"):
        if smi_col is None:
            smi_col = self.smi_col
        if act_col is None:
            act_col = self.act_col

        data = self.data
        smi_list = data[smi_col].tolist()
        act_list = data[act_col].tolist()
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

    def get_cliffs(self):
        return pd.DataFrame(self.act_list)  # a df containing smiles string, sali score, ic50 and molecule


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
