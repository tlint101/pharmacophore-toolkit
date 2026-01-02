set_color Donor_color, (0.2549019607843137, 0.4117647058823529, 0.8823529411764706)
set_color Acceptor_color, (1.0, 0.27058823529411763, 0.0)
set_color Aromatic_color, (0.8549019607843137, 0.6470588235294118, 0.12549019607843137)
set_color Hydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color LumpedHydrophobe_color, (0.1803921568627451, 0.5450980392156862, 0.3411764705882353)
set_color PosIonizable_color, (0.0, 0.7490196078431373, 1.0)
pseudoatom Donor_1, pos=[2.9960301450011877, -0.25162248165329726, 1.1343162349988476]
pseudoatom Acceptor_1, pos=[0.0003826897260829215, -4.193965755119221, 0.7917336837971352]
pseudoatom Acceptor_2, pos=[-2.8992823806372265, -0.494938213006445, 0.13590536838610062]
pseudoatom Acceptor_3, pos=[-2.5282940117940567, -3.057594735270693, 0.9310719642162649]
pseudoatom Aromatic_1, pos=[-0.37109267257338985, -1.5859987793957364, -0.018812249206457455]
pseudoatom Hydrophobe_1, pos=[-3.3169031057808906, 0.7859784721700365, -0.21303076308279423]
pseudoatom Hydrophobe_2, pos=[-2.9415094573274487, -3.125333696952897, 2.2756687417290378]
pseudoatom Hydrophobe_3, pos=[1.3145400100499454, -4.755033081101116, 0.7091775062029062]

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
