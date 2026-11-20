"""
Hold constant values
"""
import os
import re
import matplotlib.colors as mcolors
from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig
from typing import Optional, Union

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


def sanitize_pymol_name(name: Optional[str], fallback: str = "Mol") -> str:
    """
    Sanitize a molecule or feature name for safe use in PyMOL object identifiers and filenames.
    :param name: Optional[str]
        Raw name string.
    :param fallback: str
        Fallback identifier returned when sanitization yields an empty string. Defaults to ``'Mol'``.
    :return: str
        Identifier containing only ASCII letters, digits, and underscores.
    """
    text = "" if name is None else str(name).strip()
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_")
    return cleaned or fallback


def unique_pymol_names(names: list[str]) -> list[str]:
    """
    Sanitize a list of molecule names and disambiguate case-insensitive duplicates so PyMOL
    objects and output filenames never collide.
    :param names: list[str]
        Sequence of raw molecule names.
    :return: list[str]
        Unique sanitized identifiers preserving the input order.
    """
    sanitized = [sanitize_pymol_name(n, fallback=f"Mol_{i}") for i, n in enumerate(names)]
    counts: dict[str, int] = {}
    for s in sanitized:
        key = s.lower()
        counts[key] = counts.get(key, 0) + 1

    used_lower: set[str] = set()
    unique: list[str] = []
    for i, s in enumerate(sanitized):
        candidate = s
        if counts[s.lower()] > 1 or candidate.lower() in used_lower:
            candidate = f"{s}_{i}"
            k = 1
            while candidate.lower() in used_lower:
                candidate = f"{s}_{i}_{k}"
                k += 1
        used_lower.add(candidate.lower())
        unique.append(candidate)
    return unique


def resolve_pml_output_paths(
    savepath: Optional[Union[str, os.PathLike]],
    mol_names: list[str],
    suffix: str,
) -> tuple[list[str], list[str]]:
    """
    Resolve one or more PyMOL ``.pml`` output paths from ``savepath`` and molecule names.

    For a single molecule, ``savepath`` may be either an explicit ``.pml`` file path or an output
    directory (writing ``<dir>/<mol_name>_<suffix>.pml``). For multiple molecules, ``savepath`` must
    be an output directory and produces one ``<mol_name>_<suffix>.pml`` file per molecule.
    :param savepath: Optional[Union[str, os.PathLike]]
        Required ``.pml`` file path (single molecule) or output directory path.
    :param mol_names: list[str]
        Molecule names corresponding to the items being exported.
    :param suffix: str
        Filename suffix placed before ``.pml`` when writing inside a directory (e.g. ``'pli'``,
        ``'cliffs'``, or ``'pharmacophore'``).
    :return: tuple[list[str], list[str]]
        Tuple of ``(output_paths, safe_names)`` aligned with ``mol_names``.
    """
    if savepath is None or not isinstance(savepath, (str, os.PathLike)):
        raise ValueError("savepath must be a non-empty .pml file path or output directory.")
    raw_path = os.fspath(savepath)
    if not isinstance(raw_path, str):
        raise ValueError("savepath must be a non-empty .pml file path or output directory.")
    savepath_str = raw_path.strip()
    if not savepath_str:
        raise ValueError("savepath must be a non-empty .pml file path or output directory.")
    if not mol_names:
        raise ValueError("At least one molecule is required to export a .pml script.")

    safe_names = unique_pymol_names(mol_names)
    is_pml_file = savepath_str.lower().endswith(".pml")

    if len(safe_names) == 1:
        if is_pml_file:
            if os.path.isdir(savepath_str):
                raise ValueError(
                    f"savepath {savepath_str!r} is an existing directory; expected a .pml file path."
                )
            parent_dir = os.path.dirname(savepath_str)
            if parent_dir:
                os.makedirs(parent_dir, exist_ok=True)
            return [savepath_str], safe_names
        if os.path.isfile(savepath_str):
            raise ValueError(
                f"savepath {savepath_str!r} is an existing file; expected a .pml file path or output directory."
            )
        os.makedirs(savepath_str, exist_ok=True)
        return [os.path.join(savepath_str, f"{safe_names[0]}_{suffix}.pml")], safe_names

    if is_pml_file:
        raise ValueError(
            f"Multiple molecules ({len(safe_names)}) require an output directory for savepath, "
            f"got .pml file path {savepath_str!r}."
        )
    if os.path.isfile(savepath_str):
        raise ValueError(
            f"savepath {savepath_str!r} is an existing file; expected an output directory for "
            f"{len(safe_names)} molecules."
        )
    os.makedirs(savepath_str, exist_ok=True)
    return [os.path.join(savepath_str, f"{name}_{suffix}.pml") for name in safe_names], safe_names


def is_hbond_acceptor_nitrogen(atom: Chem.Atom) -> bool:
    """
    Return True when a nitrogen atom can act as a hydrogen-bond acceptor (works with both explicit
    and implicit hydrogens). Excludes positively charged/4-coordinate N, pyrrole/indole-like
    aromatic [nX3], conjugated arylamines ([NX3]-a), and neutral amide/sulfonamide/phosphoramide/
    amidine nitrogens while preserving deprotonated anions (e.g. sulfonamide [N-]), pyridines,
    imines, nitriles, and neutral aliphatic amines.
    :param atom: Chem.Atom
        RDKit atom to inspect.
    :return: bool
        True if the atom is a valid hydrogen-bond acceptor nitrogen.
    """
    if atom.GetAtomicNum() != 7:
        return False
    chg = atom.GetFormalCharge()
    if chg > 0:
        return False
    if chg < 0:
        return True
    total_deg = atom.GetTotalDegree()
    if total_deg >= 4:
        return False
    if atom.GetIsAromatic():
        return total_deg < 3
    if total_deg == 3:
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() != 1.0:
                continue
            nb = bond.GetOtherAtom(atom)
            if nb.GetIsAromatic():
                return False
            if nb.GetAtomicNum() in (6, 15, 16) and any(
                b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(nb).GetAtomicNum() in (7, 8, 16)
                for b in nb.GetBonds()
            ):
                return False
    return True


def format_pymol_split_spheres_block(
    split_specs: list[tuple[str, float, float, float, tuple[float, float, float], float, float]],
    sphere_scale: Optional[float] = None,
    transparency: float = 0.0,
) -> str:
    """
    Generate a PyMOL python...python end CGO block rendering multi-color spherical
    wedges for coincident pharmacophore features so PyMOL matches 3D View split spheres.
    :param split_specs: list
        Entries of (obj_name, x, y, z, (r, g, b), phi0, phi1).
    :param sphere_scale: Optional[float]
        Multiplier applied to PyMOL's pseudoatom base VdW radius (None uses 1.0x base VdW).
    :param transparency: float
        CGO transparency in [0.0, 1.0] (0.0 opaque).
    :return: str
        PyMOL python...python end script block, or '' when split_specs is empty.
    """
    if not split_specs:
        return ""

    scale_expr = f"_ph_base_r * {float(sphere_scale)}" if sphere_scale is not None else "_ph_base_r"
    transp_val = round(float(transparency), 3)
    lines = [
        "python",
        "import math as _math",
        "from pymol import cmd as _cmd, cgo as _cgo",
        "_ph_prev_az = _cmd.get('auto_zoom')",
        "_cmd.set('auto_zoom', 0, quiet=1)",
        "_cmd.pseudoatom('_ph_vdw_probe', pos=[0.0, 0.0, 0.0], quiet=1)",
        "_ph_mdl = _cmd.get_model('_ph_vdw_probe')",
        "_ph_base_r = float(_ph_mdl.atom[0].vdw) if (_ph_mdl.atom and float(_ph_mdl.atom[0].vdw) > 0) else 1.0",
        "_cmd.delete('_ph_vdw_probe')",
        "_cmd.set('auto_zoom', _ph_prev_az, quiet=1)",
        f"_ph_r = {scale_expr}",
        "def _ph_load_wedge(name, cx, cy, cz, r, rgb, phi0, phi1, transp):",
        "    n_lat, n_lon = 20, 16",
        "    stride = n_lon + 1",
        "    verts, norms = [], []",
        "    for i in range(n_lat + 1):",
        "        th = _math.pi * i / n_lat",
        "        st, ct = _math.sin(th), _math.cos(th)",
        "        for j in range(n_lon + 1):",
        "            ph = phi0 + (phi1 - phi0) * j / n_lon",
        "            nx, ny, nz = st * _math.sin(ph), ct, st * _math.cos(ph)",
        "            verts.append((cx + r * nx, cy + r * ny, cz + r * nz))",
        "            norms.append((nx, ny, nz))",
        "    obj = [_cgo.BEGIN, _cgo.TRIANGLES, _cgo.COLOR, rgb[0], rgb[1], rgb[2]]",
        "    for i in range(n_lat):",
        "        for j in range(n_lon):",
        "            v00 = i * stride + j",
        "            v10 = (i + 1) * stride + j",
        "            v11 = (i + 1) * stride + (j + 1)",
        "            v01 = i * stride + (j + 1)",
        "            if i > 0:",
        "                for k in (v00, v11, v01):",
        "                    obj.extend([_cgo.NORMAL, *norms[k], _cgo.VERTEX, *verts[k]])",
        "            if i < n_lat - 1:",
        "                for k in (v00, v10, v11):",
        "                    obj.extend([_cgo.NORMAL, *norms[k], _cgo.VERTEX, *verts[k]])",
        "    obj.append(_cgo.END)",
        "    _cmd.load_cgo(obj, name)",
        "    if transp > 0.0:",
        "        _cmd.set('cgo_transparency', transp, name)",
    ]
    for obj_name, x, y, z, rgb, phi0, phi1 in split_specs:
        r_c, g_c, b_c = float(rgb[0]), float(rgb[1]), float(rgb[2])
        lines.append(
            f"_ph_load_wedge({obj_name!r}, {x:.5f}, {y:.5f}, {z:.5f}, _ph_r, "
            f"({r_c:.4f}, {g_c:.4f}, {b_c:.4f}), {phi0:.6f}, {phi1:.6f}, {transp_val})"
        )
    lines.append("python end")
    return "\n".join(lines) + "\n"


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

# Default 3D interaction cutoffs (Å and degrees, Maestro-aligned)
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
DEFAULT_PI_FACE_DIST = 4.4
DEFAULT_PI_FACE_ANGLE = 30.0
DEFAULT_PI_EDGE_ANGLE = 60.0
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
    "HBond_Donor": (0.2549019607843137, 0.4117647058823529, 0.8823529411764706),     # royalblue (when split_hbonds=True)
    "HBond_Acceptor": (1.0, 0.27058823529411763, 0.0),                              # orangered (when split_hbonds=True)
    "Hydrophobic": (0.39215686274509803, 0.4549019607843137, 0.5450980392156862),   # slate gray (#64748B)
    "Pi_Stacking": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),   # goldenrod
    "Cation_Pi": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),     # goldenrod
    "Ionic": (0.7294117647058823, 0.3333333333333333, 0.8274509803921568),          # mediumorchid
    "Clash": (0.8627450980392157, 0.0784313725490196, 0.23529411764705882),         # crimson
}

PML_COLORS = {
    "HBond": (0.0, 1.0, 0.0),                                                          # PyMOL native green
    "HBond_Donor": (0.2549019607843137, 0.4117647058823529, 0.8823529411764706),       # royalblue (when split_hbonds=True)
    "HBond_Acceptor": (1.0, 0.27058823529411763, 0.0),                                # orangered (when split_hbonds=True)
    "Hydrophobic": (0.1803921568627451, 0.5450980392156862, 0.3411764705882353),      # seagreen
    "Pi_Stacking": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),     # goldenrod
    "Cation_Pi": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),       # goldenrod
    "Ionic": (0.7294117647058823, 0.3333333333333333, 0.8274509803921568),            # mediumorchid
    "Clash": (0.8627450980392157, 0.0784313725490196, 0.23529411764705882),           # crimson
}

PROLIG_COLORS = {
    "HBDonor": (0.2549019607843137, 0.4117647058823529, 0.8823529411764706),        # royalblue (#4169E1)
    "HBAcceptor": (1.0, 0.27058823529411763, 0.0),                                   # orangered (#FF4500)
    "XBDonor": (0.050980392156862744, 0.5803921568627451, 0.5333333333333333),       # teal (#0D9488)
    "XBAcceptor": (0.023529411764705882, 0.7137254901960784, 0.8313725490196079),    # cyan (#06B6D4)
    "PiStacking": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),     # goldenrod (#DAA520)
    "FaceToFace": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),     # goldenrod (#DAA520)
    "EdgeToFace": (0.7058823529411765, 0.3254901960784314, 0.03529411764705882),     # dark amber (#B45309)
    "CationPi": (0.8509803921568627, 0.4666666666666667, 0.023529411764705882),      # amber (#D97706)
    "PiCation": (0.9607843137254902, 0.6196078431372549, 0.043137254901960784),      # gold-orange (#F59E0B)
    "Cationic": (0.0, 0.7490196078431373, 1.0),                                      # deepskyblue (#00BFFF)
    "Anionic": (0.7294117647058823, 0.3333333333333333, 0.8274509803921568),         # mediumorchid (#BA55D3)
    "MetalDonor": (0.27450980392156865, 0.5098039215686274, 0.7058823529411765),     # steelblue (#4682B4)
    "MetalAcceptor": (0.38823529411764707, 0.4, 0.9450980392156862),                 # indigo (#6366F1)
    "Hydrophobic": (0.39215686274509803, 0.4549019607843137, 0.5450980392156862),    # slate gray (#64748B)
    "VdWContact": (0.796078431372549, 0.8352941176470589, 0.8823529411764706),       # light cool gray (#CBD5E1)
}


if __name__ == "__main__":
    import doctest

    doctest.testmod()
