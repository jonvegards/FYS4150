e_var = [4.44071, 4.41601, 4.60681, 4.80239, 5.03183, 5.40081, 5.50442, 5.80788, 6.08185, 6.45361, 7.02942, 7.33624, 7.69243, 8.56567, 8.98515, 10.1885, 11.4666, 12.4718, 12.5706, 11.9861, 10.8073, 9.81194, 9.10693, 8.19381, 7.91946, 7.63205, 7.32048, 6.98433, 6.77534, 6.6245, 6.48845, 6.33738, 6.15385, 6.09307, 5.96519, 5.84133, 5.724, 5.64038, 5.55259, 5.46475, ];
m_var = [6.68685, 2.22476, 2.49746, 2.93788, 3.29132, 5.2333, 4.44651, 5.20656, 6.538, 8.4506, 14.1256, 14.2817, 16.9297, 34.087, 38.2753, 81.551, 157.578, 216.061, 296.399, 345.808, 320.802, 278.697, 227.866, 170.896, 151.966, 126.995, 100.316, 79.8638, 73.4892, 62.2318, 55.5243, 51.5489, 42.801, 39.8522, 36.2713, 33.0053, 30.1275, 27.5239, 26.2461, 23.7731, ];
mean_mag = [0.86167, 0.856366, 0.850586, 0.843716, 0.836914, 0.834434, 0.826136, 0.817702, 0.807567, 0.796329, 0.798717, 0.784985, 0.769283, 0.747912, 0.727407, 0.674173, 0.624719, 0.567134, 0.490887, 0.405312, 0.333952, 0.276523, 0.231844, 0.187697, 0.175448, 0.158138, 0.139344, 0.123226, 0.117584, 0.107152, 0.101054, 0.0963852, 0.0880155, 0.0843701, 0.0805319, 0.0769025, 0.0732396, 0.0700424, 0.0685223, 0.0651394, ];
temperatures = [2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, ];
test = e_var./temperatures;
plot(temperatures,test)
