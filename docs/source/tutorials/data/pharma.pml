set_color Donor_color, (0.2549019607843137, 0.4117647058823529, 0.8823529411764706)
set_color Acceptor_color, (1.0, 0.27058823529411763, 0.0)
set_color Aromatic_color, (0.8549019607843137, 0.6470588235294118, 0.12549019607843137)
set_color Hydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color LumpedHydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color PosIonizable_color, (0.0, 0.7490196078431373, 1.0)
pseudoatom Donor_1, pos=[2.897378451604482, -0.33046028896054835, 1.1758990001045737]
pseudoatom Acceptor_1, pos=[-1.3425815336054603, -2.990171846612501, 0.8116986128481385]
pseudoatom Acceptor_2, pos=[-2.570876703483008, 1.406277539492896, -0.4720299733117135]
pseudoatom Acceptor_3, pos=[-3.197678369972157, -1.0339967739196383, 0.55368084698354]
pseudoatom Aromatic_1, pos=[-0.6172431522975114, -0.5394469281673627, -0.24150057548229553]
pseudoatom Hydrophobe_1, pos=[-2.4589141550504214, 2.70300457491222, -0.9640723842684452]
pseudoatom Hydrophobe_2, pos=[-4.231690974271513, -1.560946169199195, -0.21033963680920187]
pseudoatom Hydrophobe_3, pos=[-0.3719057935567369, -4.02309348049242, 0.9497118598179433]

show spheres, Acceptor_*
color acceptor_color, Acceptor_*

show spheres, Donor_*
color donor_color, Donor_*

show spheres, Hydrophobe_*
color hydrophobe_color, Hydrophobe_*

show spheres, Aromatic_*
color aromatic_color, Aromatic_*

show spheres, LumpedHydrophobe_*
color lumpedhydrophobe, LumpedHydrophobe_*

show spheres, PosIonizable*
color posionizable, PosIonizable*

set sphere_scale, 0.7
