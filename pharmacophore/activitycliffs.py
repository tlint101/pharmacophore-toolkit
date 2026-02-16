import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from typing import Optional


class ActivityCliffs:  #todo finalize class name
    def __init__(self, smi: list, activity: list):
        self.smi_list = smi
        self.act_list = activity

    def get_sali(self, ic50: bool = False):
        if ic50 is True:
            pass  # convert ic50 to pic50 level

    def get_cliffs(self):
        return pd.DataFrame(self.act_list)  # a df containing smiles string, sali score, ic50 and molecule

    def calc_LiPE(self, smi: Optional[list], activity: Optional[list], units: str = "nM"):  #todo test code
        """
        LipE = pIC50 - LogP
        https://en.wikipedia.org/wiki/Lipophilic_efficiency
        :return:
        """
        if smi is None:
            smi_list = self.smi_list
        if activity is None:
            act_list = self.act_list

        # calculate ic50
        pIC = _calculate_pic50(act_list, units)

        # calculate LogP
        logp = []
        for x in smi_list:
            mol = Chem.MolFromSmiles(x)
            calc_logp = MolLogP(mol)
            logp.append(calc_logp)

        # calculate LiPE
        LiPE = []
        for x, y in zip(pIC, logp):
            calc_LiPE = x - y
            LiPE.append(calc_LiPE)

        return LiPE


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

        pic50.append(converted_ic)

    return pic50
