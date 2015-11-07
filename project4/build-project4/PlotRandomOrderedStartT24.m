% Comparing numerical results, ordered and random start state

MC = [100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1e+06, 1.1e+06, 1.2e+06, 1.3e+06, 1.4e+06, 1.5e+06, 1.6e+06, 1.7e+06, 1.8e+06, 1.9e+06, 2e+06, 2.1e+06, 2.2e+06, 2.3e+06, 2.4e+06, 2.5e+06, 2.6e+06, 2.7e+06, 2.8e+06, 2.9e+06, 3e+06, 3.1e+06, 3.2e+06, 3.3e+06, 3.4e+06, 3.5e+06, 3.6e+06, 3.7e+06, 3.8e+06, 3.9e+06, 4e+06, 4.1e+06, 4.2e+06, 4.3e+06, 4.4e+06, 4.5e+06, 4.6e+06, 4.7e+06, 4.8e+06, 4.9e+06, 5e+06, 5.1e+06, 5.2e+06, 5.3e+06, 5.4e+06, 5.5e+06, 5.6e+06, 5.7e+06, 5.8e+06, 5.9e+06, 6e+06, 6.1e+06, 6.2e+06, 6.3e+06, 6.4e+06, 6.5e+06, 6.6e+06, 6.7e+06, 6.8e+06, 6.9e+06, 7e+06, 7.1e+06, 7.2e+06, 7.3e+06, 7.4e+06, 7.5e+06, 7.6e+06, 7.7e+06, 7.8e+06, 7.9e+06, 8e+06, 8.1e+06, 8.2e+06, 8.3e+06, 8.4e+06, 8.5e+06, 8.6e+06, 8.7e+06, 8.8e+06, 8.9e+06, 9e+06, 9.1e+06, 9.2e+06, 9.3e+06, 9.4e+06, 9.5e+06, 9.6e+06, 9.7e+06, 9.8e+06, 9.9e+06, 1e+07, ];
e_var = [0.415297, 0.413213, 0.413597, 0.412946, 0.413472, 0.41301, 0.413122, 0.412832, 0.413136, 0.413215, 0.413433, 0.413315, 0.413001, 0.412909, 0.412702, 0.412639, 0.41242, 0.412718, 0.412666, 0.412428, 0.412442, 0.41255, 0.412585, 0.412521, 0.412389, 0.412107, 0.412375, 0.412243, 0.412169, 0.412254, 0.412357, 0.412422, 0.412454, 0.41239, 0.412403, 0.412406, 0.412508, 0.412558, 0.412517, 0.412575, 0.412616, 0.412613, 0.412547, 0.412479, 0.412484, 0.412504, 0.412546, 0.412583, 0.412588, 0.412659, 0.412714, 0.412666, 0.412663, 0.412704, 0.412776, 0.412791, 0.412929, 0.412777, 0.412725, 0.412712, 0.412698, 0.412844, 0.412927, 0.412927, 0.412907, 0.412998, 0.412898, 0.412819, 0.412775, 0.412756, 0.412746, 0.412697, 0.412693, 0.412576, 0.412546, 0.412605, 0.412573, 0.412627, 0.412655, 0.412684, 0.412687, 0.412762, 0.412789, 0.412793, 0.412835, 0.412839, 0.412874, 0.412849, 0.41285, 0.412825, 0.412806, 0.412741, 0.412745, 0.412723, 0.412623, 0.412648, 0.412666, 0.412663, 0.412589, 0.412604, ];
m_var = [0.125938, 0.125082, 0.125464, 0.125267, 0.125589, 0.125453, 0.125588, 0.125594, 0.125427, 0.125326, 0.125434, 0.125293, 0.125274, 0.125223, 0.125171, 0.12514, 0.124953, 0.125055, 0.124951, 0.124883, 0.124884, 0.124903, 0.124964, 0.124914, 0.12488, 0.124797, 0.124912, 0.124924, 0.12486, 0.12486, 0.124896, 0.124956, 0.124994, 0.124945, 0.124946, 0.124966, 0.124986, 0.12498, 0.124987, 0.125012, 0.125004, 0.124985, 0.124968, 0.124964, 0.124972, 0.125013, 0.125051, 0.125063, 0.125081, 0.125094, 0.125099, 0.125088, 0.125073, 0.125108, 0.125136, 0.125151, 0.125187, 0.125129, 0.125115, 0.125091, 0.125106, 0.125145, 0.125172, 0.125177, 0.125174, 0.125202, 0.125159, 0.12514, 0.125111, 0.125109, 0.125099, 0.125089, 0.125068, 0.12503, 0.125037, 0.125058, 0.125059, 0.125078, 0.125091, 0.125099, 0.125096, 0.125123, 0.12513, 0.125124, 0.125137, 0.125142, 0.125153, 0.125139, 0.125155, 0.125146, 0.125143, 0.125131, 0.125139, 0.125139, 0.125112, 0.125111, 0.12512, 0.125125, 0.125101, 0.125106, ];
mean_mag = [0.88068, 0.881428, 0.881137, 0.88133, 0.88103, 0.881215, 0.881106, 0.881183, 0.881218, 0.881213, 0.881138, 0.88124, 0.881327, 0.881375, 0.88145, 0.881488, 0.881626, 0.881513, 0.881584, 0.881676, 0.881684, 0.881675, 0.881643, 0.881685, 0.881733, 0.881836, 0.881723, 0.88174, 0.881783, 0.881763, 0.881732, 0.88169, 0.881658, 0.881689, 0.881685, 0.88168, 0.881655, 0.881653, 0.881659, 0.881636, 0.881635, 0.881651, 0.881673, 0.881694, 0.881686, 0.881665, 0.881641, 0.881629, 0.881622, 0.881594, 0.881578, 0.881584, 0.881585, 0.881557, 0.881529, 0.881518, 0.881473, 0.881529, 0.881543, 0.88156, 0.881557, 0.881512, 0.881481, 0.881477, 0.881487, 0.881455, 0.881495, 0.881518, 0.88154, 0.881545, 0.881552, 0.881563, 0.881574, 0.881615, 0.881618, 0.881597, 0.881605, 0.88158, 0.881564, 0.881551, 0.881549, 0.881523, 0.881515, 0.88152, 0.881508, 0.881504, 0.881491, 0.881506, 0.881497, 0.881506, 0.88151, 0.881528, 0.881525, 0.881531, 0.881562, 0.88156, 0.881551, 0.881549, 0.881574, 0.881569, ];
time_used = [0.066963, 0.134778, 0.20055, 0.268123, 0.333466, 0.399639, 0.466456, 0.532221, 0.597722, 0.665371, 0.729495, 0.795952, 0.862064, 0.927067, 0.992564, 1.06082, 1.12812, 1.19397, 1.25958, 1.32517, 1.38961, 1.45505, 1.51939, 1.58484, 1.65152, 1.71792, 1.7823, 1.85148, 1.91941, 1.98432, 2.04953, 2.11577, 2.18158, 2.24799, 2.31453, 2.38037, 2.44609, 2.51133, 2.57699, 2.64149, 2.70727, 2.77303, 2.84054, 2.90378, 2.969, 3.03455, 3.09967, 3.16507, 3.23029, 3.29651, 3.36228, 3.42926, 3.49452, 3.55981, 3.62426, 3.68964, 3.75463, 3.82205, 3.88922, 3.95565, 4.02077, 4.08738, 4.15405, 4.22028, 4.28661, 4.35435, 4.4199, 4.48659, 4.55217, 4.61763, 4.68377, 4.74886, 4.81657, 4.88083, 4.9483, 5.01526, 5.08496, 5.15227, 5.21799, 5.28288, 5.34836, 5.4147, 5.47946, 5.5449, 5.61101, 5.67729, 5.74506, 5.81269, 5.87846, 5.94374, 6.00892, 6.07319, 6.13898, 6.20367, 6.27026, 6.33587, 6.40201, 6.46939, 6.53243, 6.59748, ];
accepted_states = [0.71616, 0.35412, 0.2396, 0.17701, 0.143716, 0.11842, 0.102564, 0.0888638, 0.0792222, 0.071798, 0.0653073, 0.059495, 0.0542, 0.05055, 0.0474013, 0.0443313, 0.0414247, 0.03995, 0.0371379, 0.035254, 0.0338429, 0.0325782, 0.0311478, 0.029475, 0.0281552, 0.0268077, 0.0268467, 0.0253143, 0.0243324, 0.0239667, 0.0231606, 0.0224819, 0.0217364, 0.0208571, 0.0203517, 0.0197114, 0.0194076, 0.0187363, 0.0182054, 0.0179288, 0.0173776, 0.0168352, 0.0165419, 0.0160814, 0.0158916, 0.0154913, 0.0152336, 0.0148771, 0.0145473, 0.0144328, 0.0140427, 0.01365, 0.0134653, 0.0133459, 0.0130553, 0.0127739, 0.0127475, 0.0120472, 0.0120358, 0.0118367, 0.0116738, 0.011724, 0.0114837, 0.0111472, 0.0108852, 0.0109794, 0.010506, 0.0104591, 0.0102446, 0.0101849, 0.01004, 0.00983833, 0.00975178, 0.00944676, 0.0094236, 0.00947118, 0.00922286, 0.00928077, 0.00911696, 0.00898625, 0.00881383, 0.00879366, 0.00866988, 0.00848786, 0.00843459, 0.00834256, 0.00826713, 0.00797795, 0.00801775, 0.00784922, 0.00784396, 0.00765935, 0.00765688, 0.00751809, 0.00732063, 0.00746938, 0.00740186, 0.00727653, 0.00706687, 0.0071914, ];
probability = [0, -8, -8, -8, 0, -8, 0, -8, -8, -8, -8, -8, -8, -8, -8, 0, -8, 0, -8, 0, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, 0, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, 0, 0, -8, -8, -8, 0, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, -8, -8, -8, -8, 0, 0, -8, 0, -8, -8, -8, 0, 0, 0, -8, -8, -8, -8, -8, -8, -8, -8, ];
avg_energy = [-1.64104, -1.64291, -1.64239, -1.64295, -1.64224, -1.64283, -1.64263, -1.64301, -1.64276, -1.64253, -1.64238, -1.64255, -1.64293, -1.64305, -1.6433, -1.6434, -1.64365, -1.64332, -1.64342, -1.64371, -1.64374, -1.64371, -1.64369, -1.64378, -1.64394, -1.64426, -1.64396, -1.64408, -1.64415, -1.64405, -1.64397, -1.6439, -1.64384, -1.64388, -1.64387, -1.64388, -1.6438, -1.64377, -1.64381, -1.64375, -1.64372, -1.64375, -1.64382, -1.64392, -1.6439, -1.64388, -1.64385, -1.64382, -1.64382, -1.64372, -1.64365, -1.64367, -1.64364, -1.64359, -1.64352, -1.6435, -1.64335, -1.6435, -1.64354, -1.64357, -1.64359, -1.64345, -1.64335, -1.64335, -1.64339, -1.64329, -1.64339, -1.64347, -1.64351, -1.64353, -1.64354, -1.64358, -1.64358, -1.6437, -1.64374, -1.64368, -1.64372, -1.64364, -1.6436, -1.64356, -1.64354, -1.64347, -1.64345, -1.64346, -1.64342, -1.64341, -1.64337, -1.64341, -1.64341, -1.64343, -1.64345, -1.64351, -1.64351, -1.64354, -1.64364, -1.64362, -1.6436, -1.6436, -1.64367, -1.64366, ];

ran_e_var = [0.415244, 0.413186, 0.413579, 0.412932, 0.413461, 0.413004, 0.413116, 0.412825, 0.413132, 0.41321, 0.413428, 0.413311, 0.412997, 0.412905, 0.412699, 0.412637, 0.412417, 0.412716, 0.412663, 0.412426, 0.412439, 0.412547, 0.412582, 0.412519, 0.412387, 0.412105, 0.412374, 0.412241, 0.412167, 0.412252, 0.412355, 0.41242, 0.412453, 0.412388, 0.412401, 0.412405, 0.412507, 0.412557, 0.412515, 0.412574, 0.412615, 0.412612, 0.412545, 0.412478, 0.412483, 0.412503, 0.412545, 0.412582, 0.412587, 0.412658, 0.412713, 0.412665, 0.412662, 0.412703, 0.412775, 0.41279, 0.412928, 0.412777, 0.412724, 0.412711, 0.412698, 0.412843, 0.412927, 0.412926, 0.412906, 0.412997, 0.412898, 0.412819, 0.412774, 0.412755, 0.412746, 0.412697, 0.412692, 0.412575, 0.412545, 0.412604, 0.412573, 0.412626, 0.412654, 0.412683, 0.412687, 0.412761, 0.412788, 0.412792, 0.412834, 0.412838, 0.412874, 0.412848, 0.41285, 0.412825, 0.412806, 0.41274, 0.412745, 0.412722, 0.412623, 0.412648, 0.412665, 0.412663, 0.412589, 0.412604, ];
ran_m_var = [0.125921, 0.125074, 0.125458, 0.125262, 0.125585, 0.125452, 0.125585, 0.125592, 0.125425, 0.125324, 0.125433, 0.125292, 0.125272, 0.125221, 0.12517, 0.125139, 0.124952, 0.125055, 0.12495, 0.124883, 0.124883, 0.124902, 0.124963, 0.124914, 0.124879, 0.124797, 0.124911, 0.124923, 0.124859, 0.124859, 0.124896, 0.124956, 0.124994, 0.124944, 0.124946, 0.124965, 0.124986, 0.12498, 0.124986, 0.125011, 0.125003, 0.124985, 0.124967, 0.124964, 0.124972, 0.125012, 0.12505, 0.125063, 0.12508, 0.125094, 0.125099, 0.125088, 0.125072, 0.125108, 0.125136, 0.125151, 0.125187, 0.125129, 0.125115, 0.125091, 0.125106, 0.125145, 0.125171, 0.125177, 0.125174, 0.125202, 0.125159, 0.12514, 0.125111, 0.125109, 0.125099, 0.125089, 0.125068, 0.12503, 0.125037, 0.125058, 0.125059, 0.125078, 0.125091, 0.125099, 0.125096, 0.125123, 0.12513, 0.125124, 0.125137, 0.125141, 0.125153, 0.125139, 0.125155, 0.125146, 0.125143, 0.125131, 0.125139, 0.125139, 0.125111, 0.125111, 0.12512, 0.125125, 0.125101, 0.125106, ];
ran_mean_mag = [0.8807, 0.881437, 0.881143, 0.881335, 0.881034, 0.881217, 0.881109, 0.881186, 0.88122, 0.881216, 0.88114, 0.881241, 0.881328, 0.881376, 0.881451, 0.881488, 0.881627, 0.881513, 0.881585, 0.881677, 0.881685, 0.881676, 0.881644, 0.881686, 0.881734, 0.881837, 0.881724, 0.881741, 0.881784, 0.881764, 0.881732, 0.88169, 0.881659, 0.88169, 0.881686, 0.881681, 0.881655, 0.881654, 0.88166, 0.881637, 0.881636, 0.881651, 0.881674, 0.881695, 0.881686, 0.881665, 0.881641, 0.881629, 0.881622, 0.881594, 0.881578, 0.881585, 0.881585, 0.881558, 0.88153, 0.881518, 0.881473, 0.881529, 0.881543, 0.881561, 0.881557, 0.881512, 0.881481, 0.881478, 0.881487, 0.881456, 0.881496, 0.881519, 0.88154, 0.881545, 0.881552, 0.881564, 0.881574, 0.881615, 0.881618, 0.881597, 0.881605, 0.88158, 0.881564, 0.881551, 0.88155, 0.881523, 0.881516, 0.881521, 0.881508, 0.881504, 0.881491, 0.881506, 0.881497, 0.881506, 0.88151, 0.881528, 0.881525, 0.881531, 0.881563, 0.88156, 0.881551, 0.881549, 0.881574, 0.881569, ];
ran_accepted_states = [0.71607, 0.35412, 0.239607, 0.177005, 0.143724, 0.118417, 0.102567, 0.0888588, 0.0792256, 0.071795, 0.0653073, 0.059495, 0.0542, 0.05055, 0.0474013, 0.0443325, 0.0414235, 0.0399522, 0.0371358, 0.035255, 0.0338419, 0.0325782, 0.0311478, 0.029475, 0.0281552, 0.0268077, 0.026847, 0.0253139, 0.0243324, 0.0239667, 0.023161, 0.0224816, 0.0217367, 0.0208568, 0.020352, 0.0197111, 0.0194076, 0.0187363, 0.0182062, 0.017928, 0.0173778, 0.016835, 0.0165423, 0.0160809, 0.0158916, 0.0154913, 0.0152336, 0.0148771, 0.0145473, 0.0144336, 0.014042, 0.01365, 0.0134653, 0.0133459, 0.0130553, 0.0127739, 0.0127482, 0.0120471, 0.0120353, 0.0118367, 0.0116743, 0.0117237, 0.0114835, 0.0111472, 0.0108852, 0.0109794, 0.0105061, 0.010459, 0.0102446, 0.0101849, 0.0100401, 0.00983819, 0.00975178, 0.00944676, 0.009424, 0.00947092, 0.00922273, 0.00928077, 0.00911696, 0.00898625, 0.00881383, 0.00879366, 0.00867012, 0.0084881, 0.00843412, 0.00834291, 0.00826678, 0.00797795, 0.00801775, 0.00784933, 0.00784429, 0.00765935, 0.00765645, 0.00751809, 0.00732074, 0.00746927, 0.00740186, 0.00727673, 0.00706667, 0.0071914, ];
ran_probability = [0, -8, -8, -8, -8, 0, 0, -8, 0, -8, -8, -8, -8, -8, -8, 0, -8, 0, -8, 0, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, 0, -8, 0, -8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, -8, -8, -8, 0, 0, -8, -8, 0, -8, -8, -8, -8, -8, 0, -8, -8, -8, 0, -8, -8, -8, -8, 0, -8, -8, -8, -8, -8, -8, -8, -8, -8, 0, -8, -8, -8, -8, 0, 0, -8, -8, 0, -8, -8, 0, -8, -8, ];
ran_avg_energy = [-1.6411, -1.64294, -1.64241, -1.64296, -1.64225, -1.64283, -1.64263, -1.64302, -1.64276, -1.64254, -1.64239, -1.64255, -1.64293, -1.64305, -1.6433, -1.6434, -1.64365, -1.64332, -1.64342, -1.64372, -1.64374, -1.64371, -1.64369, -1.64378, -1.64394, -1.64426, -1.64396, -1.64408, -1.64415, -1.64405, -1.64397, -1.6439, -1.64384, -1.64389, -1.64387, -1.64389, -1.6438, -1.64377, -1.64381, -1.64375, -1.64373, -1.64375, -1.64382, -1.64392, -1.6439, -1.64388, -1.64385, -1.64382, -1.64382, -1.64372, -1.64365, -1.64367, -1.64364, -1.64359, -1.64352, -1.6435, -1.64335, -1.6435, -1.64355, -1.64357, -1.64359, -1.64345, -1.64335, -1.64335, -1.64339, -1.64329, -1.6434, -1.64347, -1.64351, -1.64353, -1.64354, -1.64358, -1.64359, -1.64371, -1.64374, -1.64368, -1.64372, -1.64364, -1.6436, -1.64356, -1.64354, -1.64347, -1.64345, -1.64346, -1.64342, -1.64341, -1.64337, -1.64341, -1.64341, -1.64343, -1.64345, -1.64351, -1.64351, -1.64354, -1.64364, -1.64362, -1.6436, -1.6436, -1.64367, -1.64366, ];

std(e_var)
std(m_var)
std(mean_mag)
std(avg_energy)

%Plotting results
figure(1)
plot(MC, m_var, MC, ran_m_var)
legend('Ordered','Random')
title('Susceptibilitet \chi som funksjon av MC-sykluser, T=2.4','FontSize',16)
xlabel('Ant. Monte Carlo-sykluser','FontSize',16)
ylabel('Susceptibilitet','FontSize',16)
set(0,'DefaultAxesFontSize',15)

figure(2)
plot(MC, e_var, MC, ran_e_var)
legend('Ordered','Random')
title('Varmekapasitet C_V som funksjon av MC-sykluser, T=2.4','FontSize',16)
xlabel('Ant. Monte Carlo-sykluser','FontSize',16)
ylabel('Varmekapasitet','FontSize',16)

figure(3)
plot(MC,avg_energy, MC, ran_avg_energy)
legend('Ordered','Random')
title('Gjennomsnittlig energi som funksjon av MC-sykluser, T=2.4','FontSize',16)
xlabel('Ant. Monte Carlo-sykluser','FontSize',16)
ylabel('\langle E \rangle','FontSize',16)
