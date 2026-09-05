"""
Hold constant values
"""
import os
import matplotlib.colors as mcolors
from rdkit.Chem import AllChem, RDConfig
from typing import Optional

feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

FEATURES = {
    "Donor": ["[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",  # nitrogen in rings or amines
              "[#8!H0&!$([OH][C,S,P]=O)]",  # oxygen atom bonded to hydrogen
              "[#16!H0]"  # sulfur atom
              ],
    "Acceptor": ["[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
                 # amines or nitrogen in rings non aromatic rings
                 "[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]"  # hydroxyl or ether
                 ],
    "Aromatic": ["a1aaaaa1",  # six membered aromatic ring
                 "a1aaaa1",  # five membered aromatic ring
                 "[#6]1[#6]=[#6][#6]=[#7]1",  # fused pyrole ring
                 "[#6]:1:[#6]:[#6]:[#6]:[#6]:[#7]:1"  # fused pyridine ring
                 ],
    "Hydrophobe": ["a1aaaaa1",  # six member aromatic ring
                   "a1aaaa1",  # five member aromatic ring
                   "*~1~*~*~*~*~*~1",  # fix member ring, any atom, any bond
                   "*~1~*~*~*~*~1",  # five member ring, any atom, any bond
                   # methyl, methylene, methine, halogens
                   "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
                   # matches to methyl, methylene, terminal methine
                   "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
                   "[CH3]",  # terminal methyl group
                   "[CH2]~*~!@[*1]",  # terminal methylene group not in a ring
                   ]
}

FEATURE_COLORS = {
    "Donor": (0.2549019607843137, 0.4117647058823529, 0.8823529411764706),  # royalblue
    "Acceptor": (1.0, 0.27058823529411763, 0.0),  # orangered
    "Aromatic": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),  # goldenrod
    "Hydrophobe": (0.1803921568627451, 0.5450980392156862, 0.3411764705882353),  # seagreen
    "LumpedHydrophobe": (0.1803921568627451, 0.5450980392156862, 0.3411764705882353),  # seagreen
    "PosIonizable": (0.0, 0.7490196078431373, 1.0)  # deepskyblue
}

# for rendering py3dmol
INTERACTIVE_COLORS = {
    "Donor": "royalblue",
    "Acceptor": "orangered",
    "Aromatic": "goldenrod",
    "Hydrophobe": "seagreen",
    "LumpedHydrophobe": "seagreen",
    "PosIonizable": "deepskyblue"
}


def color_convert(color: Optional[str] = None):
    """
    helper function to convert color to rgb.
    :param color: Optional[str]
        Color to convert to rgb. Can be hex or color name.
    :return:
    """
    try:
        # convert color to rgb
        rgb = mcolors.to_rgb(color)
        return rgb
    except:
        raise ValueError(f"{color} is not a valid color!")


if __name__ == "__main__":
    import doctest

    doctest.testmod()
