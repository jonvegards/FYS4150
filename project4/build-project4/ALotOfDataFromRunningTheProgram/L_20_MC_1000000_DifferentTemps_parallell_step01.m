e_var = [7.35676, 7.55931, 7.95055, 8.12845, 8.4292, 8.69666, 9.06711, 9.23317, 9.32717, 9.52446, ];
m_var = [7.97811, 8.64725, 9.52586, 10.3329, 11.3709, 12.7036, 14.1492, 15.3807, 16.0679, 17.4041, ];
mean_mag = [0.785023, 0.787981, 0.776001, 0.745086, 0.749104, 0.729277, 0.706021, 0.694557, 0.676068, 0.659974, ];
temperatures = [2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, ];
test = e_var./temperatures;
plot(temperatures,m_var)