%two electron case_armadillo.m
%% Two electron case, several states
R1= [0.00286718, 0.00580006, 0.00879406, 0.0118445, 0.0149464, 0.018095, 0.021285, 0.0245114, 0.0277688, 0.031052, 0.0343555, 0.037674, 0.0410019, 0.0443338, 0.0476642, 0.0509876, 0.0542986, 0.0575918, 0.0608618, 0.0641034, 0.0673112, 0.0704801, 0.0736052, 0.0766815, 0.0797041, 0.0826685, 0.0855701, 0.0884045, 0.0911675, 0.0938551, 0.0964635, 0.0989889, 0.101428, 0.103777, 0.106034, 0.108195, 0.110258, 0.11222, 0.114079, 0.115834, 0.117482, 0.119022, 0.120453, 0.121773, 0.122983, 0.12408, 0.125066, 0.125939, 0.1267, 0.12735, 0.127887, 0.128315, 0.128632, 0.128841, 0.128943, 0.128939, 0.128831, 0.12862, 0.128309, 0.1279, 0.127395, 0.126796, 0.126106, 0.125328, 0.124464, 0.123517, 0.122491, 0.121388, 0.120211, 0.118964, 0.11765, 0.116271, 0.114832, 0.113336, 0.111785, 0.110184, 0.108536, 0.106844, 0.105111, 0.103341, 0.101537, 0.0997029, 0.0978411, 0.0959552, 0.0940481, 0.0921231, 0.0901831, 0.088231, 0.0862698, 0.0843022, 0.082331, 0.0803588, 0.0783881, 0.0764214, 0.074461, 0.0725093, 0.0705684, 0.0686404, 0.0667272, 0.0648308, 0.0629529, 0.0610952, 0.0592594, 0.0574468, 0.0556589, 0.0538969, 0.0521621, 0.0504556, 0.0487783, 0.0471312, 0.0455151, 0.0439307, 0.0423787, 0.0408596, 0.0393739, 0.037922, 0.0365043, 0.035121, 0.0337722, 0.0324582, 0.0311789, 0.0299343, 0.0287244, 0.027549, 0.026408, 0.0253011, 0.024228, 0.0231885, 0.0221821, 0.0212085, 0.0202672, 0.0193577, 0.0184796, 0.0176323, 0.0168152, 0.0160279, 0.0152696, 0.0145399, 0.0138379, 0.0131632, 0.012515, 0.0118927, 0.0112956, 0.0107231, 0.0101744, 0.00964892, 0.00914593, 0.00866476, 0.00820472, 0.00776515, 0.00734538, 0.00694474, 0.00656258, 0.00619825, 0.00585111, 0.00552052, 0.00520589, 0.00490658, 0.00462202, 0.00435161, 0.00409477, 0.00385096, 0.00361963, 0.00340025, 0.00319229, 0.00299525, 0.00280864, 0.00263199, 0.00246483, 0.00230671, 0.0021572, 0.00201587, 0.00188232, 0.00175615, 0.00163699, 0.00152445, 0.00141819, 0.00131786, 0.00122312, 0.00113367, 0.00104918, 0.000969349, 0.0008939, 0.000822545, 0.000755014, 0.000691052, 0.000630395, 0.000572798, 0.000518026, 0.000465842, 0.000416018, 0.000368331, 0.000322564, 0.000278504, 0.000235926, 0.000194636, 0.000154412, 0.000115049, 7.63326e-05, 3.80587e-05, ];
R2= [0.00396578, 0.00801585, 0.0121371, 0.016316, 0.0205389, 0.0247916, 0.0290601, 0.0333299, 0.0375866, 0.0418157, 0.0460028, 0.0501335, 0.0541936, 0.0581689, 0.0620457, 0.0658104, 0.0694498, 0.0729511, 0.0763018, 0.0794902, 0.0825048, 0.0853348, 0.0879699, 0.0904006, 0.0926181, 0.094614, 0.096381, 0.0979123, 0.0992019, 0.100245, 0.101036, 0.101574, 0.101853, 0.101874, 0.101634, 0.101135, 0.100375, 0.099357, 0.098083, 0.0965561, 0.0947802, 0.0927597, 0.0905, 0.0880071, 0.0852878, 0.0823495, 0.0792, 0.0758479, 0.0723023, 0.0685729, 0.0646696, 0.060603, 0.056384, 0.0520238, 0.0475338, 0.0429259, 0.038212, 0.0334044, 0.0285152, 0.0235568, 0.0185416, 0.013482, 0.00839038, 0.00327892, -0.00184021, -0.00695507, -0.0120539, -0.0171251, -0.0221576, -0.0271402, -0.0320625, -0.036914, -0.041685, -0.046366, -0.0509479, -0.0554222, -0.0597806, -0.0640157, -0.0681201, -0.0720874, -0.0759113, -0.0795862, -0.0831072, -0.0864696, -0.0896694, -0.0927032, -0.0955681, -0.0982614, -0.100781, -0.103127, -0.105296, -0.10729, -0.109107, -0.110749, -0.112216, -0.113509, -0.114632, -0.115584, -0.11637, -0.116992, -0.117453, -0.117757, -0.117907, -0.117908, -0.117763, -0.117478, -0.117057, -0.116506, -0.115828, -0.11503, -0.114116, -0.113092, -0.111964, -0.110737, -0.109417, -0.108008, -0.106518, -0.104952, -0.103314, -0.101611, -0.0998482, -0.0980307, -0.0961639, -0.0942531, -0.0923033, -0.0903195, -0.0883065, -0.0862691, -0.0842119, -0.0821391, -0.0800552, -0.0779641, -0.0758697, -0.073776, -0.0716863, -0.0696042, -0.0675329, -0.0654753, -0.0634345, -0.061413, -0.0594134, -0.057438, -0.0554891, -0.0535686, -0.0516783, -0.04982, -0.0479951, -0.046205, -0.044451, -0.042734, -0.041055, -0.0394148, -0.037814, -0.0362532, -0.0347327, -0.0332528, -0.0318138, -0.0304156, -0.0290583, -0.0277418, -0.0264657, -0.0252299, -0.0240339, -0.0228773, -0.0217597, -0.0206803, -0.0196386, -0.0186339, -0.0176655, -0.0167326, -0.0158343, -0.0149698, -0.0141381, -0.0133385, -0.0125699, -0.0118313, -0.0111218, -0.0104403, -0.00978585, -0.00915734, -0.00855374, -0.00797397, -0.00741696, -0.00688163, -0.00636687, -0.00587158, -0.00539465, -0.00493498, -0.00449142, -0.00406285, -0.00364811, -0.00324607, -0.00285555, -0.00247536, -0.00210431, -0.00174119, -0.00138475, -0.00103373, -0.000686823, -0.000342701, ];
R3= [0.00479105, 0.00967575, 0.0146298, 0.0196285, 0.0246464, 0.0296582, 0.0346381, 0.0395604, 0.0443996, 0.0491305, 0.053728, 0.058168, 0.0624266, 0.066481, 0.0703093, 0.0738904, 0.0772047, 0.0802337, 0.08296, 0.0853682, 0.0874438, 0.0891745, 0.0905494, 0.0915592, 0.0921966, 0.092456, 0.0923338, 0.0918282, 0.090939, 0.0896684, 0.08802, 0.0859995, 0.0836143, 0.0808735, 0.077788, 0.0743704, 0.0706349, 0.0665969, 0.0622736, 0.0576834, 0.0528458, 0.0477817, 0.0425128, 0.0370617, 0.031452, 0.0257077, 0.0198537, 0.013915, 0.007917, 0.00188537, -0.00415427, -0.0101763, -0.0161554, -0.0220664, -0.0278847, -0.0335862, -0.0391475, -0.0445458, -0.0497593, -0.0547671, -0.0595493, -0.0640872, -0.0683633, -0.0723612, -0.0760659, -0.0794639, -0.0825429, -0.0852922, -0.0877024, -0.0897658, -0.0914762, -0.0928287, -0.09382, -0.0944486, -0.094714, -0.0946175, -0.0941618, -0.0933509, -0.0921901, -0.0906862, -0.0888471, -0.086682, -0.084201, -0.0814157, -0.0783383, -0.074982, -0.071361, -0.0674901, -0.0633848, -0.0590615, -0.0545366, -0.0498273, -0.0449512, -0.039926, -0.0347696, -0.0295002, -0.0241357, -0.0186944, -0.0131942, -0.00765284, -0.00208789, 0.00348342, 0.00904424, 0.0145782, 0.0200692, 0.0255021, 0.0308618, 0.0361344, 0.0413062, 0.0463644, 0.051297, 0.0560925, 0.0607404, 0.0652309, 0.0695551, 0.0737048, 0.0776725, 0.0814518, 0.085037, 0.0884232, 0.0916061, 0.0945826, 0.09735, 0.0999065, 0.102251, 0.104384, 0.106304, 0.108013, 0.109513, 0.110806, 0.111895, 0.112782, 0.113473, 0.11397, 0.114279, 0.114405, 0.114353, 0.114129, 0.11374, 0.113191, 0.112489, 0.111641, 0.110653, 0.109533, 0.108288, 0.106925, 0.105452, 0.103874, 0.1022, 0.100437, 0.0985919, 0.0966713, 0.0946822, 0.0926316, 0.0905258, 0.0883714, 0.0861745, 0.0839412, 0.0816773, 0.0793883, 0.0770798, 0.0747568, 0.0724242, 0.0700867, 0.0677486, 0.0654141, 0.0630871, 0.0607713, 0.0584698, 0.056186, 0.0539225, 0.0516821, 0.0494669, 0.0472792, 0.0451206, 0.0429929, 0.0408973, 0.038835, 0.0368068, 0.0348134, 0.0328553, 0.0309326, 0.0290453, 0.0271934, 0.0253763, 0.0235934, 0.0218441, 0.0201273, 0.0184418, 0.0167864, 0.0151596, 0.0135595, 0.0119844, 0.0104323, 0.00890083, 0.00738769, 0.00589029, 0.00440586, 0.00293147, 0.00146396, ];
R4= [0.00576613, 0.0116343, 0.0175643, 0.0235152, 0.0294457, 0.0353141, 0.0410792, 0.0466998, 0.0521358, 0.0573478, 0.062298, 0.0669498, 0.0712689, 0.0752227, 0.0787811, 0.0819166, 0.0846044, 0.0868226, 0.0885524, 0.0897783, 0.0904882, 0.0906733, 0.0903285, 0.089452, 0.0880459, 0.0861156, 0.0836702, 0.0807221, 0.0772873, 0.073385, 0.0690375, 0.0642701, 0.059111, 0.0535909, 0.0477431, 0.0416029, 0.0352077, 0.0285965, 0.0218098, 0.0148891, 0.00787699, 0.000816504, -0.00624899, -0.0132761, -0.0202218, -0.0270435, -0.0336996, -0.0401495, -0.0463539, -0.0522753, -0.057878, -0.0631285, -0.0679952, -0.0724495, -0.0764651, -0.0800186, -0.0830896, -0.0856607, -0.0877176, -0.0892492, -0.0902478, -0.0907089, -0.0906311, -0.0900165, -0.0888703, -0.087201, -0.08502, -0.0823417, -0.0791835, -0.0755653, -0.0715098, -0.0670419, -0.062189, -0.0569803, -0.0514468, -0.0456214, -0.0395382, -0.0332326, -0.0267408, -0.0200998, -0.0133474, -0.00652115, 0.000340901, 0.00720108, 0.0140221, 0.0207672, 0.0274004, 0.033887, 0.040193, 0.0462862, 0.0521356, 0.0577122, 0.0629885, 0.067939, 0.0725404, 0.0767714, 0.080613, 0.0840482, 0.0870628, 0.0896445, 0.0917838, 0.0934732, 0.0947078, 0.0954851, 0.0958047, 0.0956687, 0.0950813, 0.0940489, 0.0925799, 0.0906846, 0.0883755, 0.0856664, 0.082573, 0.0791126, 0.0753038, 0.0711664, 0.0667214, 0.0619909, 0.0569977, 0.0517654, 0.0463181, 0.0406806, 0.0348778, 0.0289347, 0.0228766, 0.0167286, 0.0105157, 0.00426235, -0.00200715, -0.00826913, -0.0145005, -0.020679, -0.0267829, -0.0327917, -0.0386855, -0.0444456, -0.0500542, -0.0554948, -0.0607518, -0.0658111, -0.0706594, -0.0752848, -0.0796768, -0.0838257, -0.0877235, -0.091363, -0.0947384, -0.0978451, -0.10068, -0.103239, -0.105523, -0.107531, -0.109263, -0.110722, -0.111909, -0.112829, -0.113485, -0.113883, -0.114029, -0.113928, -0.113587, -0.113015, -0.112219, -0.111207, -0.109988, -0.108571, -0.106965, -0.105179, -0.103224, -0.101108, -0.0988416, -0.0964339, -0.0938947, -0.0912334, -0.0884594, -0.0855817, -0.0826094, -0.0795511, -0.0764152, -0.0732099, -0.069943, -0.0666218, -0.0632535, -0.0598448, -0.0564018, -0.0529304, -0.049436, -0.0459236, -0.0423977, -0.0388622, -0.0353206, -0.0317761, -0.0282312, -0.024688, -0.021148, -0.0176123, -0.0140814, -0.0105555, -0.00703398, -0.00351599, ];
R5= [0.00703298, 0.0141749, 0.021361, 0.0285255, 0.0356024, 0.042526, 0.0492314, 0.0556554, 0.0617368, 0.0674171, 0.0726412, 0.0773575, 0.081519, 0.0850833, 0.0880131, 0.0902767, 0.0918481, 0.0927076, 0.0928415, 0.0922428, 0.0909108, 0.0888515, 0.0860772, 0.0826069, 0.0784656, 0.0736845, 0.0683004, 0.0623557, 0.0558977, 0.0489786, 0.0416544, 0.0339853, 0.0260342, 0.0178668, 0.00955086, 0.00115538, -0.00724973, -0.0155944, -0.023809, -0.031825, -0.0395755, -0.0469959, -0.0540242, -0.060602, -0.0666746, -0.0721917, -0.0771076, -0.081382, -0.0849798, -0.0878717, -0.0900345, -0.091451, -0.0921105, -0.0920084, -0.0911466, -0.0895335, -0.0871833, -0.0841167, -0.08036, -0.0759453, -0.0709097, -0.0652957, -0.0591502, -0.052524, -0.0454721, -0.0380524, -0.0303257, -0.0223548, -0.0142044, -0.00594022, 0.0023714, 0.0106641, 0.0188721, 0.0269306, 0.0347763, 0.0423482, 0.0495875, 0.0564386, 0.0628494, 0.0687713, 0.0741603, 0.0789763, 0.0831843, 0.0867543, 0.0896612, 0.0918853, 0.0934125, 0.0942337, 0.0943456, 0.0937502, 0.0924547, 0.0904718, 0.0878191, 0.0845188, 0.0805982, 0.0760886, 0.0710256, 0.0654482, 0.0593992, 0.0529241, 0.0460713, 0.0388912, 0.0314364, 0.0237605, 0.0159185, 0.00796566, -4.22437e-05, -0.00804987, -0.0160024, -0.023846, -0.0315283, -0.0389986, -0.046208, -0.0531104, -0.0596619, -0.0658217, -0.0715522, -0.076819, -0.0815912, -0.0858414, -0.0895464, -0.0926863, -0.0952455, -0.0972122, -0.0985785, -0.0993403, -0.0994976, -0.099054, -0.0980166, -0.0963964, -0.0942075, -0.0914673, -0.0881961, -0.0844174, -0.0801569, -0.0754428, -0.0703056, -0.0647777, -0.0588929, -0.0526867, -0.0461957, -0.0394573, -0.0325096, -0.0253913, -0.0181409, -0.0107972, -0.00339859, 0.00401719, 0.0114129, 0.0187524, 0.0260005, 0.0331234, 0.0400884, 0.0468648, 0.0534233, 0.0597364, 0.0657785, 0.0715261, 0.0769576, 0.0820534, 0.0867962, 0.0911707, 0.0951638, 0.0987646, 0.101964, 0.104755, 0.107134, 0.109097, 0.110643, 0.111774, 0.112493, 0.112803, 0.112711, 0.112224, 0.111352, 0.110104, 0.108492, 0.106528, 0.104226, 0.101599, 0.0986638, 0.0954346, 0.0919278, 0.0881599, 0.0841475, 0.0799075, 0.0754569, 0.0708126, 0.0659916, 0.0610104, 0.0558855, 0.0506329, 0.0452682, 0.0398065, 0.0342624, 0.0286498, 0.0229821, 0.0172718, 0.0115308, 0.00577009, ];
rho= [0, 0.0247525, 0.049505, 0.0742574, 0.0990099, 0.123762, 0.148515, 0.173267, 0.19802, 0.222772, 0.247525, 0.272277, 0.29703, 0.321782, 0.346535, 0.371287, 0.39604, 0.420792, 0.445545, 0.470297, 0.49505, 0.519802, 0.544554, 0.569307, 0.594059, 0.618812, 0.643564, 0.668317, 0.693069, 0.717822, 0.742574, 0.767327, 0.792079, 0.816832, 0.841584, 0.866337, 0.891089, 0.915842, 0.940594, 0.965347, 0.990099, 1.01485, 1.0396, 1.06436, 1.08911, 1.11386, 1.13861, 1.16337, 1.18812, 1.21287, 1.23762, 1.26238, 1.28713, 1.31188, 1.33663, 1.36139, 1.38614, 1.41089, 1.43564, 1.4604, 1.48515, 1.5099, 1.53465, 1.55941, 1.58416, 1.60891, 1.63366, 1.65842, 1.68317, 1.70792, 1.73267, 1.75743, 1.78218, 1.80693, 1.83168, 1.85644, 1.88119, 1.90594, 1.93069, 1.95545, 1.9802, 2.00495, 2.0297, 2.05446, 2.07921, 2.10396, 2.12871, 2.15347, 2.17822, 2.20297, 2.22772, 2.25248, 2.27723, 2.30198, 2.32673, 2.35149, 2.37624, 2.40099, 2.42574, 2.4505, 2.47525, 2.5, 2.52475, 2.5495, 2.57426, 2.59901, 2.62376, 2.64851, 2.67327, 2.69802, 2.72277, 2.74752, 2.77228, 2.79703, 2.82178, 2.84653, 2.87129, 2.89604, 2.92079, 2.94554, 2.9703, 2.99505, 3.0198, 3.04455, 3.06931, 3.09406, 3.11881, 3.14356, 3.16832, 3.19307, 3.21782, 3.24257, 3.26733, 3.29208, 3.31683, 3.34158, 3.36634, 3.39109, 3.41584, 3.44059, 3.46535, 3.4901, 3.51485, 3.5396, 3.56436, 3.58911, 3.61386, 3.63861, 3.66337, 3.68812, 3.71287, 3.73762, 3.76238, 3.78713, 3.81188, 3.83663, 3.86139, 3.88614, 3.91089, 3.93564, 3.9604, 3.98515, 4.0099, 4.03465, 4.05941, 4.08416, 4.10891, 4.13366, 4.15842, 4.18317, 4.20792, 4.23267, 4.25743, 4.28218, 4.30693, 4.33168, 4.35644, 4.38119, 4.40594, 4.43069, 4.45545, 4.4802, 4.50495, 4.5297, 4.55446, 4.57921, 4.60396, 4.62871, 4.65347, 4.67822, 4.70297, 4.72772, 4.75248, 4.77723, 4.80198, 4.82673, 4.85149, 4.87624, 4.90099, 4.92574, 4.9505, 4.97525, ];
% figure(1)
% plot(rho(2:201),R1.*R1)
% hold('on')
% plot(rho(2:201),R2.*R2)
% plot(rho(2:201),R3.*R3)
% plot(rho(2:201),R4.*R4)
% plot(rho(2:201),R5.*R5)
% legend('R1','R2','R3','R4','R5')
% xlabel('Distance \rho')
% ylabel('Probability |\Psi|^2')

%% Two electron case
R1_8= [-0.00584322, -0.0118904, -0.018117, -0.0244978, -0.0310064, -0.0376155, -0.0442976, -0.0510244, -0.0577678, -0.0644992, -0.0711904, -0.0778135, -0.0843409, -0.0907458, -0.0970023, -0.103085, -0.10897, -0.114636, -0.120059, -0.125221, -0.130103, -0.134689, -0.138963, -0.142913, -0.146528, -0.149798, -0.152715, -0.155275, -0.157473, -0.159308, -0.16078, -0.161891, -0.162644, -0.163045, -0.1631, -0.162817, -0.162207, -0.161281, -0.16005, -0.158527, -0.156728, -0.154666, -0.152358, -0.149819, -0.147067, -0.144119, -0.140991, -0.137702, -0.134269, -0.130708, -0.127038, -0.123276, -0.119437, -0.115537, -0.111593, -0.10762, -0.103631, -0.0996408, -0.0956618, -0.0917064, -0.087786, -0.0839113, -0.080092, -0.0763371, -0.0726547, -0.0690523, -0.0655363, -0.0621124, -0.0587857, -0.0555604, -0.0524398, -0.0494268, -0.0465236, -0.0437316, -0.0410518, -0.0384845, -0.0360295, -0.0336863, -0.0314537, -0.0293303, -0.0273142, -0.0254034, -0.0235955, -0.0218876, -0.0202771, -0.0187608, -0.0173355, -0.0159978, -0.0147445, -0.013572, -0.0124769, -0.0114555, -0.0105044, -0.00962009, -0.00879912, -0.00803809, -0.00733367, -0.00668259, -0.00608172, -0.00552797, -0.00501839, -0.00455012, -0.00412044, -0.00372673, -0.00336647, -0.00303729, -0.00273693, -0.00246325, -0.00221421, -0.00198792, -0.00178258, -0.0015965, -0.0014281, -0.00127591, -0.00113856, -0.00101476, -0.000903327, -0.00080316, -0.000713238, -0.000632621, -0.00056044, -0.000495899, -0.000438263, -0.000386861, -0.000341079, -0.000300356, -0.000264178, -0.000232081, -0.000203641, -0.000178473, -0.000156229, -0.000136595, -0.000119287, -0.000104048, -9.06491e-05, -7.88818e-05, -6.85609e-05, -5.95199e-05, -5.16102e-05, -4.46989e-05, -3.86674e-05, -3.34105e-05, -2.88343e-05, -2.48556e-05, -2.14008e-05, -1.84045e-05, -1.58092e-05, -1.3564e-05, -1.1624e-05, -9.94984e-06, -8.50684e-06, -7.26462e-06, -6.19655e-06, -5.27935e-06, -4.49267e-06, -3.81876e-06, -3.24216e-06, -2.74942e-06, -2.32886e-06, -1.97033e-06, -1.66507e-06, -1.40547e-06, -1.18496e-06, -9.97895e-07, -8.39387e-07, -7.0524e-07, -5.91845e-07, -4.96109e-07, -4.15379e-07, -3.47383e-07, -2.90183e-07, -2.42121e-07, -2.01787e-07, -1.67977e-07, -1.39671e-07, -1.16e-07, -9.62296e-08, -7.9736e-08, -6.59925e-08, -5.45538e-08, -4.50446e-08, -3.71484e-08, -3.05991e-08, -2.51729e-08, -2.06819e-08, -1.69686e-08, -1.39013e-08, -1.13694e-08, -9.28062e-09, -7.55792e-09, -6.13684e-09, -4.96353e-09, -3.99295e-09, -3.18735e-09, -2.51496e-09, -1.94887e-09, -1.46609e-09, -1.0467e-09, -6.73058e-10, -3.29173e-10, ];
R2_8= [0.00807036, 0.0163881, 0.0248837, 0.0334852, 0.0421191, 0.0507108, 0.0591854, 0.0674685, 0.0754868, 0.0831688, 0.0904458, 0.097252, 0.103526, 0.109209, 0.11425, 0.118601, 0.122222, 0.125076, 0.127135, 0.128377, 0.128787, 0.128357, 0.127085, 0.124976, 0.122043, 0.118305, 0.113786, 0.108517, 0.102534, 0.0958795, 0.0885991, 0.0807434, 0.0723663, 0.063525, 0.0542793, 0.0446907, 0.0348224, 0.0247382, 0.0145022, 0.00417827, -0.00617075, -0.0164832, -0.0266995, -0.0367619, -0.0466158, -0.0562095, -0.0654946, -0.0744263, -0.082964, -0.091071, -0.0987146, -0.105867, -0.112504, -0.118606, -0.124159, -0.129152, -0.133578, -0.137434, -0.140722, -0.143447, -0.145618, -0.147246, -0.148345, -0.148934, -0.149032, -0.148661, -0.147844, -0.146608, -0.144978, -0.142982, -0.140648, -0.138006, -0.135083, -0.13191, -0.128516, -0.124927, -0.121172, -0.117279, -0.113273, -0.109178, -0.10502, -0.10082, -0.0966002, -0.09238, -0.088178, -0.0840111, -0.079895, -0.0758437, -0.0718699, -0.0679848, -0.0641982, -0.0605187, -0.0569534, -0.0535083, -0.0501883, -0.0469969, -0.0439368, -0.0410098, -0.0382166, -0.0355573, -0.033031, -0.0306364, -0.0283714, -0.0262336, -0.0242198, -0.0223268, -0.0205508, -0.0188878, -0.0173336, -0.0158838, -0.0145339, -0.0132793, -0.0121154, -0.0110376, -0.0100412, -0.00912169, -0.00827461, -0.00749558, -0.00678032, -0.00612472, -0.00552478, -0.00497667, -0.00447674, -0.00402147, -0.00360754, -0.00323178, -0.00289122, -0.00258302, -0.00230456, -0.00205334, -0.00182704, -0.00162351, -0.00144072, -0.0012768, -0.00113004, -0.000998814, -0.000881664, -0.000777229, -0.000684265, -0.000601631, -0.000528284, -0.000463275, -0.000405738, -0.000354886, -0.000310005, -0.000270451, -0.00023564, -0.000205045, -0.000178194, -0.000154661, -0.000134065, -0.000116063, -0.000100351, -8.66551e-05, -7.47341e-05, -6.43714e-05, -5.53757e-05, -4.75771e-05, -4.08254e-05, -3.49878e-05, -2.99474e-05, -2.5601e-05, -2.18581e-05, -1.86391e-05, -1.58744e-05, -1.35029e-05, -1.14715e-05, -9.73357e-06, -8.24873e-06, -6.98175e-06, -5.90207e-06, -4.9832e-06, -4.20219e-06, -3.53922e-06, -2.97717e-06, -2.50129e-06, -2.09888e-06, -1.75903e-06, -1.47237e-06, -1.23089e-06, -1.02771e-06, -8.5697e-07, -7.13658e-07, -5.93509e-07, -4.92889e-07, -4.08712e-07, -3.38355e-07, -2.79594e-07, -2.30543e-07, -1.89607e-07, -1.55432e-07, -1.26872e-07, -1.02955e-07, -8.28543e-08, -6.58646e-08, -5.13792e-08, -3.8872e-08, -2.78801e-08, -1.79887e-08, -8.81607e-09, ];
R3_8= [-0.00958525, -0.0194229, -0.0293875, -0.0393504, -0.0491814, -0.0587508, -0.0679309, -0.0765979, -0.0846335, -0.0919268, -0.0983755, -0.103887, -0.108382, -0.111791, -0.114059, -0.115147, -0.115027, -0.11369, -0.111139, -0.107395, -0.102491, -0.096475, -0.0894102, -0.0813709, -0.0724436, -0.0627254, -0.0523225, -0.0413488, -0.0299247, -0.0181753, -0.00622879, 0.00578505, 0.0177366, 0.0294977, 0.0409437, 0.0519544, 0.0624155, 0.0722202, 0.08127, 0.0894758, 0.0967588, 0.103051, 0.108297, 0.112451, 0.115482, 0.11737, 0.118106, 0.117696, 0.116154, 0.113508, 0.109794, 0.10506, 0.0993607, 0.0927604, 0.0853299, 0.0771458, 0.06829, 0.0588483, 0.0489093, 0.0385635, 0.0279021, 0.0170162, 0.00599582, -0.00507109, -0.0160994, -0.0270075, -0.0377181, -0.0481587, -0.0582621, -0.067967, -0.077218, -0.085966, -0.0941686, -0.10179, -0.1088, -0.115177, -0.120904, -0.12597, -0.13037, -0.134107, -0.137186, -0.139619, -0.141421, -0.142612, -0.143216, -0.143258, -0.142768, -0.141778, -0.14032, -0.13843, -0.136143, -0.133495, -0.130523, -0.127262, -0.12375, -0.120021, -0.116108, -0.112047, -0.107867, -0.103599, -0.0992712, -0.0949109, -0.0905424, -0.0861886, -0.0818706, -0.0776071, -0.0734153, -0.0693102, -0.0653052, -0.0614116, -0.0576391, -0.0539957, -0.050488, -0.047121, -0.0438982, -0.0408221, -0.0378938, -0.0351135, -0.0324804, -0.0299928, -0.0276483, -0.0254438, -0.0233758, -0.0214401, -0.0196322, -0.0179473, -0.0163803, -0.014926, -0.0135791, -0.0123341, -0.0111855, -0.010128, -0.00915621, -0.00826484, -0.00744877, -0.00670301, -0.00602274, -0.00540334, -0.00484035, -0.00432953, -0.00386687, -0.00344853, -0.00307092, -0.00273066, -0.00242455, -0.00214963, -0.00190313, -0.00168247, -0.00148526, -0.0013093, -0.00115254, -0.00101311, -0.000889297, -0.000779516, -0.00068233, -0.000596427, -0.000520614, -0.000453808, -0.000395028, -0.000343388, -0.00029809, -0.000258412, -0.000223711, -0.000193405, -0.000166979, -0.000143968, -0.000123961, -0.000106591, -9.15313e-05, -7.84943e-05, -6.72242e-05, -5.74954e-05, -4.91091e-05, -4.18902e-05, -3.5685e-05, -3.03585e-05, -2.57926e-05, -2.18842e-05, -1.85432e-05, -1.56909e-05, -1.32592e-05, -1.11888e-05, -9.42817e-06, -7.93285e-06, -6.66431e-06, -5.58932e-06, -4.6792e-06, -3.90924e-06, -3.25817e-06, -2.70771e-06, -2.24211e-06, -1.84783e-06, -1.5132e-06, -1.22813e-06, -9.83866e-07, -7.72764e-07, -5.88079e-07, -4.2378e-07, -2.74378e-07, -1.34755e-07, ];[0.00782934, 0.0158504, 0.0239852, 0.0321534, 0.0402737, 0.0482647, 0.0560453, 0.0635358, 0.0706588, 0.0773399, 0.0835088, 0.0890994, 0.0940512, 0.0983094, 0.101826, 0.104559, 0.106476, 0.10755, 0.107764, 0.107106, 0.105577, 0.103183, 0.0999392, 0.0958683, 0.0910014, 0.0853769, 0.0790403, 0.0720435, 0.0644445, 0.0563067, 0.0476984, 0.0386917, 0.0293623, 0.0197884, 0.0100501, 0.000228522, -0.00959496, -0.0193394, -0.0289251, -0.0382743, -0.0473119, -0.055966, -0.0641688, -0.0718568, -0.0789715, -0.0854598, -0.0912745, -0.0963744, -0.100725, -0.104297, -0.10707, -0.109028, -0.110164, -0.110475, -0.109967, -0.10865, -0.106543, -0.103667, -0.100051, -0.0957284, -0.0907367, -0.085118, -0.0789178, -0.0721849, -0.0649708, -0.057329, -0.0493149, -0.0409847, -0.0323956, -0.0236048, -0.0146692, -0.00564505, 0.00341267, 0.0124504, 0.0214167, 0.0302621, 0.0389398, 0.0474058, 0.0556191, 0.0635419, 0.0711397, 0.0783816, 0.0852402, 0.0916916, 0.0977156, 0.103295, 0.108418, 0.113074, 0.117257, 0.120963, 0.124194, 0.126951, 0.129241, 0.131072, 0.132454, 0.1334, 0.133926, 0.134047, 0.133781, 0.133148, 0.132169, 0.130864, 0.129255, 0.127365, 0.125216, 0.122832, 0.120235, 0.117448, 0.114492, 0.111391, 0.108164, 0.104833, 0.101417, 0.0979359, 0.0944067, 0.0908469, 0.0872725, 0.0836986, 0.080139, 0.0766066, 0.0731133, 0.0696699, 0.0662859, 0.0629701, 0.0597302, 0.0565728, 0.0535037, 0.0505277, 0.0476489, 0.0448704, 0.0421947, 0.0396235, 0.0371578, 0.0347981, 0.0325444, 0.030396, 0.0283518, 0.0264104, 0.0245699, 0.0228282, 0.0211827, 0.0196308, 0.0181697, 0.0167962, 0.0155072, 0.0142994, 0.0131695, 0.0121141, 0.0111297, 0.010213, 0.00936056, 0.00856903, 0.00783513, 0.00715563, 0.00652738, 0.00594734, 0.00541255, 0.00492015, 0.00446739, 0.00405164, 0.00367039, 0.00332123, 0.00300187, 0.00271016, 0.00244404, 0.00220158, 0.00198094, 0.00178042, 0.0015984, 0.0014334, 0.00128398, 0.00114885, 0.00102679, 0.000916651, 0.000817389, 0.000728038, 0.000647693, 0.000575519, 0.000510759, 0.000452707, 0.000400721, 0.000354202, 0.000312613, 0.000275455, 0.000242272, 0.000212652, 0.000186219, 0.000162616, 0.000141539, 0.000122696, 0.00010582, 9.06727e-05, 7.70288e-05, 6.46868e-05, 5.34584e-05, 4.31574e-05, 3.3622e-05, 2.46909e-05, 1.62061e-05, 8.02871e-06, ];
R4_8 = [-0.0107467, -0.0217296, -0.03276, -0.043646, -0.0541958, -0.0642215, -0.0735423, -0.0819877, -0.0894006, -0.0956404, -0.100585, -0.104135, -0.106212, -0.106762, -0.105759, -0.103201, -0.0991127, -0.0935436, -0.0865691, -0.078288, -0.0688209, -0.0583082, -0.0469076, -0.0347918, -0.0221449, -0.00915989, 0.00396497, 0.0170293, 0.0298341, 0.0421848, 0.0538944, 0.064787, 0.0746996, 0.0834856, 0.0910164, 0.0971834, 0.1019, 0.105101, 0.106747, 0.106819, 0.105325, 0.102295, 0.0977805, 0.0918553, 0.0846129, 0.0761647, 0.0666381, 0.0561745, 0.0449266, 0.0330559, 0.0207306, 0.00812238, -0.00459577, -0.0172522, -0.0296792, -0.041715, -0.0532061, -0.0640093, -0.073993, -0.0830393, -0.0910449, -0.0979219, -0.103599, -0.108021, -0.111151, -0.112968, -0.113467, -0.11266, -0.110574, -0.10725, -0.102744, -0.097122, -0.090461, -0.0828479, -0.0743772, -0.0651498, -0.0552712, -0.0448501, -0.0339971, -0.022823, -0.0114377, 5.12055e-05, 0.0115397, 0.0229281, 0.0341223, 0.0450344, 0.0555833, 0.0656955, 0.0753054, 0.0843555, 0.0927967, 0.100588, 0.107698, 0.114101, 0.119782, 0.124732, 0.128949, 0.132439, 0.135213, 0.137289, 0.138689, 0.139439, 0.139571, 0.139118, 0.138117, 0.136608, 0.134631, 0.132226, 0.129438, 0.126307, 0.122876, 0.119186, 0.115276, 0.111186, 0.106953, 0.102612, 0.098196, 0.0937364, 0.0892618, 0.0847988, 0.0803716, 0.0760017, 0.0717087, 0.0675097, 0.0634193, 0.0594504, 0.0556135, 0.0519173, 0.0483686, 0.0449725, 0.0417324, 0.0386505, 0.0357275, 0.0329629, 0.0303553, 0.0279021, 0.0256001, 0.0234455, 0.0214336, 0.0195596, 0.017818, 0.0162031, 0.0147091, 0.01333, 0.0120596, 0.0108919, 0.00982076, 0.00884027, 0.00794456, 0.00712791, 0.00638481, 0.00570995, 0.00509822, 0.00454478, 0.00404499, 0.0035945, 0.00318918, 0.00282516, 0.00249881, 0.00220676, 0.00194586, 0.00171319, 0.00150605, 0.00132196, 0.00115863, 0.00101396, 0.000886036, 0.0007731, 0.000673564, 0.000585979, 0.000509035, 0.000441548, 0.00038245, 0.000330779, 0.000285674, 0.00024636, 0.000212148, 0.000182419, 0.000156627, 0.000134282, 0.000114951, 9.8253e-05, 8.38479e-05, 7.14373e-05, 6.07582e-05, 5.15793e-05, 4.36974e-05, 3.69343e-05, 3.11338e-05, 2.61593e-05, 2.18908e-05, 1.82234e-05, 1.5065e-05, 1.23344e-05, 9.96005e-06, 7.87804e-06, 6.03107e-06, 4.36694e-06, 2.83732e-06, 1.39647e-06, ];
R5_8= [0.0116953, 0.0235966, 0.0354463, 0.0469849, 0.0579564, 0.0681142, 0.0772267, 0.0850825, 0.0914952, 0.0963076, 0.0993953, 0.10067, 0.100081, 0.0976164, 0.0933051, 0.0872141, 0.0794485, 0.0701491, 0.0594894, 0.0476719, 0.0349239, 0.0214926, 0.00764025, -0.00636202, -0.0202394, -0.0337193, -0.0465371, -0.0584415, -0.0691997, -0.0786023, -0.0864671, -0.0926431, -0.0970131, -0.0994961, -0.100049, -0.0986647, -0.095377, -0.0902546, -0.0834018, -0.0749556, -0.0650828, -0.0539762, -0.0418507, -0.0289385, -0.015485, -0.00174307, 0.0120311, 0.0255836, 0.038667, 0.0510456, 0.0624991, 0.0728266, 0.0818503, 0.089418, 0.0954056, 0.0997189, 0.102294, 0.1031, 0.102136, 0.0994315, 0.0950464, 0.0890678, 0.0816086, 0.0728044, 0.0628111, 0.0518012, 0.0399608, 0.0274859, 0.0145789, 0.00144502, -0.0117111, -0.0246888, -0.0372946, -0.049345, -0.060669, -0.0711106, -0.0805308, -0.0888087, -0.0958434, -0.101554, -0.105881, -0.108787, -0.110252, -0.110279, -0.108891, -0.106126, -0.102042, -0.0967101, -0.0902165, -0.082658, -0.0741418, -0.0647826, -0.0547012, -0.044022, -0.0328718, -0.0213773, -0.00966391, 0.00214623, 0.0139354, 0.0255917, 0.0370101, 0.0480936, 0.0587536, 0.0689108, 0.0784957, 0.0874485, 0.0957196, 0.103269, 0.110067, 0.116093, 0.121336, 0.125793, 0.129469, 0.132376, 0.134535, 0.13597, 0.136712, 0.136796, 0.13626, 0.135147, 0.1335, 0.131365, 0.128789, 0.125819, 0.122502, 0.118885, 0.115012, 0.110927, 0.106673, 0.102289, 0.0978137, 0.0932815, 0.0887252, 0.0841746, 0.0796571, 0.0751968, 0.0708157, 0.0665327, 0.0623641, 0.0583241, 0.054424, 0.0506734, 0.0470794, 0.0436473, 0.0403806, 0.0372813, 0.0343498, 0.0315852, 0.0289854, 0.0265475, 0.0242675, 0.0221408, 0.0201623, 0.0183261, 0.0166262, 0.0150563, 0.0136099, 0.0122802, 0.0110606, 0.0099445, 0.00892532, 0.00799664, 0.0071522, 0.00638596, 0.00569211, 0.00506506, 0.00449952, 0.00399045, 0.00353311, 0.00312302, 0.002756, 0.00242814, 0.0021358, 0.00187561, 0.00164445, 0.00143945, 0.00125796, 0.00109757, 0.000956063, 0.000831422, 0.000721816, 0.000625581, 0.000541215, 0.000467358, 0.000402788, 0.000346402, 0.000297214, 0.000254337, 0.000216978, 0.000184427, 0.000156048, 0.000131272, 0.000109589, 9.05393e-05, 7.37087e-05, 5.87198e-05, 4.52261e-05, 3.29059e-05, 2.14556e-05, 1.05827e-05, ];
rho_8 = [0, 0.0346535, 0.0693069, 0.10396, 0.138614, 0.173267, 0.207921, 0.242574, 0.277228, 0.311881, 0.346535, 0.381188, 0.415842, 0.450495, 0.485149, 0.519802, 0.554455, 0.589109, 0.623762, 0.658416, 0.693069, 0.727723, 0.762376, 0.79703, 0.831683, 0.866337, 0.90099, 0.935644, 0.970297, 1.00495, 1.0396, 1.07426, 1.10891, 1.14356, 1.17822, 1.21287, 1.24752, 1.28218, 1.31683, 1.35149, 1.38614, 1.42079, 1.45545, 1.4901, 1.52475, 1.55941, 1.59406, 1.62871, 1.66337, 1.69802, 1.73267, 1.76733, 1.80198, 1.83663, 1.87129, 1.90594, 1.94059, 1.97525, 2.0099, 2.04455, 2.07921, 2.11386, 2.14851, 2.18317, 2.21782, 2.25248, 2.28713, 2.32178, 2.35644, 2.39109, 2.42574, 2.4604, 2.49505, 2.5297, 2.56436, 2.59901, 2.63366, 2.66832, 2.70297, 2.73762, 2.77228, 2.80693, 2.84158, 2.87624, 2.91089, 2.94554, 2.9802, 3.01485, 3.0495, 3.08416, 3.11881, 3.15347, 3.18812, 3.22277, 3.25743, 3.29208, 3.32673, 3.36139, 3.39604, 3.43069, 3.46535, 3.5, 3.53465, 3.56931, 3.60396, 3.63861, 3.67327, 3.70792, 3.74257, 3.77723, 3.81188, 3.84653, 3.88119, 3.91584, 3.9505, 3.98515, 4.0198, 4.05446, 4.08911, 4.12376, 4.15842, 4.19307, 4.22772, 4.26238, 4.29703, 4.33168, 4.36634, 4.40099, 4.43564, 4.4703, 4.50495, 4.5396, 4.57426, 4.60891, 4.64356, 4.67822, 4.71287, 4.74752, 4.78218, 4.81683, 4.85149, 4.88614, 4.92079, 4.95545, 4.9901, 5.02475, 5.05941, 5.09406, 5.12871, 5.16337, 5.19802, 5.23267, 5.26733, 5.30198, 5.33663, 5.37129, 5.40594, 5.44059, 5.47525, 5.5099, 5.54455, 5.57921, 5.61386, 5.64851, 5.68317, 5.71782, 5.75248, 5.78713, 5.82178, 5.85644, 5.89109, 5.92574, 5.9604, 5.99505, 6.0297, 6.06436, 6.09901, 6.13366, 6.16832, 6.20297, 6.23762, 6.27228, 6.30693, 6.34158, 6.37624, 6.41089, 6.44554, 6.4802, 6.51485, 6.5495, 6.58416, 6.61881, 6.65347, 6.68812, 6.72277, 6.75743, 6.79208, 6.82673, 6.86139, 6.89604, 6.93069, 6.96535, ];
figure(2)
plot(rho_8(2:201),R1_8 .*R1_8)
hold('on')
plot(rho_8(2:201),R2_8.*R2_8)
plot(rho_8(2:201),R3_8.*R3_8)
plot(rho_8(2:201),R4_8.*R4_8)
plot(rho_8(2:201),R5_8.*R5_8)
set(0,'DefaultAxesFontSize',15)
legend('\lambda_0','\lambda_1','\lambda_2','\lambda_3','\lambda_4')
title('Relative distance between two electrons for different eigenvalues, $$\omega_r = 0.5$$','FontSize',16,'Interpreter','latex')
xlabel('Distance $$\rho$$','Interpreter','latex')
ylabel('Probability $$|\Psi|^2$$','Interpreter','latex')