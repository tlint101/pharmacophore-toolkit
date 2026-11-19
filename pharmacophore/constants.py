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
    "PosIonizable": (0.0, 0.7490196078431373, 1.0),  # deepskyblue
    "NegIonizable": (0.7294117647058823, 0.3333333333333333, 0.8274509803921568),  # mediumorchid
}

# for rendering py3dmol
INTERACTIVE_COLORS = {
    "Donor": "royalblue",
    "Acceptor": "orangered",
    "Aromatic": "goldenrod",
    "Hydrophobe": "seagreen",
    "LumpedHydrophobe": "seagreen",
    "PosIonizable": "deepskyblue",
    "NegIonizable": "mediumorchid",
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


# 3-letter amino acid codes
AMINO_ACID_ALIASES = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "CYX": "CYS",
    "CYM": "CYS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "MSE": "MET",
}

STANDARD_AMINO_ACIDS = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL",
    *AMINO_ACID_ALIASES.keys(),
})

# protonation/force-field variants
WATER_RESIDUES = frozenset({"HOH", "WAT", "SOL", "DOD", "H2O"})

METAL_RESIDUES = frozenset({"MG", "ZN", "MN", "FE", "CA", "CU", "CO", "NI", "NA", "K", "CD", "HG"})

METAL_ATOMIC_NUMS = frozenset({11, 12, 19, 20, 25, 26, 27, 28, 29, 30, 48, 80})

DEFAULT_POCKET_CUTOFF = 6.0

# Default 3D interaction cutoffs (Å and degrees with relaxed offset pi-stacking)
DEFAULT_HBOND_DIST = 3.5
DEFAULT_HBOND_H_DIST = 2.8
DEFAULT_HBOND_ANGLE = 120.0
DEFAULT_HBOND_ACCEPTOR_ANGLE = 90.0
DEFAULT_HBOND_HEAVY_ANGLE = 90.0
DEFAULT_XBOND_DIST = 3.5
DEFAULT_XBOND_DONOR_ANGLE = 120.0
DEFAULT_XBOND_ACCEPTOR_ANGLE = 90.0
DEFAULT_HYDROPHOBIC_DIST = 4.0
DEFAULT_IONIC_DIST = 5.0
DEFAULT_PI_STACKING_DIST = 5.5
DEFAULT_PI_FACE_DIST = 5.0
DEFAULT_PI_FACE_ANGLE = 40.0
DEFAULT_PI_EDGE_ANGLE = 50.0
DEFAULT_CATION_PI_DIST = 6.0
DEFAULT_CATION_PI_ANGLE = 30.0
DEFAULT_METAL_DIST = 2.8
DEFAULT_VDW_TOLERANCE = 0.0
DEFAULT_CLASH_DIST = 2.2

PROLIG_INTERACTION_ORDER = (
    "HBDonor",
    "HBAcceptor",
    "XBDonor",
    "XBAcceptor",
    "PiStacking",
    "FaceToFace",
    "EdgeToFace",
    "CationPi",
    "PiCation",
    "Cationic",
    "Anionic",
    "MetalDonor",
    "MetalAcceptor",
    "Hydrophobic",
    "VdWContact",
)

PROLIG_DEFAULT_INTERACTIONS = (
    "HBDonor",
    "HBAcceptor",
    "PiStacking",
    "Hydrophobic",
    "VdWContact",
)

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

HYDROPHOBIC_RESIDUES = frozenset({
    "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR", "CYS",
})

PROTEIN_AROMATIC_RESIDUES = frozenset({"PHE", "TYR", "TRP", "HIS"})

# Standard PDB ring atom names per aromatic residue (fallback when ring perception is partial)
PROTEIN_AROMATIC_RING_ATOMS = {
    "PHE": (("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),),
    "TYR": (("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),),
    "HIS": (("CG", "ND1", "CE1", "NE2", "CD2"),),
    "TRP": (
        ("CG", "CD1", "NE1", "CE2", "CD2"),
        ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"),
    ),
}

PROTEIN_POS_ATOMS = {
    "LYS": frozenset({"NZ"}),
    "ARG": frozenset({"CZ", "NE", "NH1", "NH2"}),
}

PROTEIN_NEG_ATOMS = {
    "ASP": frozenset({"OD1", "OD2"}),
    "GLU": frozenset({"OE1", "OE2"}),
}

PLI_COLORS = {
    "HBond": (0.13333333333333333, 0.5450980392156862, 0.13333333333333333),        # forestgreen (unified default)
    "HBond_Donor": (0.13333333333333333, 0.5450980392156862, 0.13333333333333333),  # forestgreen (when split_hbonds=True)
    "HBond_Acceptor": (1.0, 0.27058823529411763, 0.0),                              # orangered (when split_hbonds=True)
    "Hydrophobic": (0.39215686274509803, 0.4549019607843137, 0.5450980392156862),   # slate gray (#64748B)
    "Pi_Stacking": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),   # goldenrod
    "Cation_Pi": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),     # goldenrod
    "Ionic": (0.7294117647058823, 0.3333333333333333, 0.8274509803921568),          # mediumorchid
    "Clash": (0.8627450980392157, 0.0784313725490196, 0.23529411764705882),         # crimson
}

PML_COLORS = {
    "HBond": (0.0, 1.0, 0.0),                                                          # PyMOL native green
    "HBond_Donor": (0.0, 1.0, 0.0),                                                    # PyMOL native green
    "HBond_Acceptor": (1.0, 0.27058823529411763, 0.0),                                # orangered (when split_hbonds=True)
    "Hydrophobic": (0.1803921568627451, 0.5450980392156862, 0.3411764705882353),      # seagreen
    "Pi_Stacking": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),     # goldenrod
    "Cation_Pi": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),       # goldenrod
    "Ionic": (0.7294117647058823, 0.3333333333333333, 0.8274509803921568),            # mediumorchid
    "Clash": (0.8627450980392157, 0.0784313725490196, 0.23529411764705882),           # crimson
}


if __name__ == "__main__":
    import doctest

    doctest.testmod()
