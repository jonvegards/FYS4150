e_var = [4.37435, 3.49492, 3.53154, 3.73768, 3.87964, 3.91185, 4.08852, 4.2228, 4.41048, 4.46316, 4.67199, 4.78382, 5.12034, 5.2105, 5.53936, 5.781, 5.84645, 6.17493, 6.87098, 7.0056, 7.25102, 8.18668, 8.14863, 8.36772, 9.09267, 9.99892, 10.3551, 10.3598, 10.6967, 10.4936, 10.775, 11.3332, 11.0332, 10.2301, 9.5222, 9.43483, 9.15534, 8.18857, 7.63501, 7.93792, 7.12633, 7.11046, 6.88647, 6.5478, 6.71015, 6.39381, 6.17519, 6.19014, 6.05915, 5.82504, 5.78374, 5.85024, 5.60345, 5.51654, 5.57481, 5.44699, 5.48638, 5.38816, 5.39216, 5.2245, 5.28289, ];
m_var = [7.51084, 1.74493, 1.75175, 2.01055, 1.94778, 2.08175, 2.30744, 2.51842, 2.51022, 2.56751, 2.83016, 3.2367, 3.59469, 3.47784, 4.28994, 5.21397, 4.34926, 5.75186, 8.61788, 10.9068, 9.02744, 18.5532, 15.469, 16.779, 29.6909, 36.9646, 37.6965, 43.8241, 49.66, 59.1304, 70.1608, 73.2185, 77.6087, 74.291, 66.1015, 70.5703, 60.6217, 54.7965, 50.0512, 52.1503, 42.2963, 39.9561, 36.3795, 35.203, 33.8637, 30.1741, 26.5029, 26.0547, 25.3249, 23.4567, 21.8839, 20.1369, 18.4449, 18.1076, 17.8191, 16.4219, 15.2818, 14.9871, 14.59, 13.9982, 14.0532, ];
mean_mag = [0.905476, 0.908369, 0.905633, 0.900562, 0.897077, 0.893246, 0.888915, 0.883499, 0.878545, 0.874829, 0.870057, 0.863551, 0.855938, 0.851236, 0.84417, 0.834765, 0.830889, 0.818791, 0.803905, 0.797929, 0.787463, 0.777155, 0.774903, 0.75115, 0.715378, 0.692622, 0.665698, 0.665045, 0.610416, 0.572872, 0.526617, 0.501616, 0.437463, 0.411869, 0.389815, 0.360138, 0.343818, 0.30723, 0.283547, 0.279783, 0.245501, 0.233513, 0.217165, 0.213031, 0.207236, 0.196264, 0.184474, 0.184302, 0.174474, 0.16963, 0.161601, 0.15362, 0.151888, 0.146712, 0.14216, 0.138858, 0.132791, 0.13027, 0.129219, 0.125778, 0.124537, ];
temperatures = [2, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.5, 2.51, 2.52, 2.53, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59, 2.6, ];
test = e_var./temperatures;
plot(temperatures,m_var./temperatures)
