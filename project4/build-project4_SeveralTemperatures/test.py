from numpy import array
import matplotlib.pylab as plt


e_var = [4.27979, 6.49511, 8.80602, 8.11945, 6.24186, ]
m_var = array([2.63431, 7.0627, 18.8094, 23.7383, 17.9669, ])
mean_mag = [0.865848, 0.792768, 0.645134, 0.460582, 0.32744, ]
temperatures = array([2.1, 2.2, 2.3, 2.4, 2.5, ])
#test = e_var./temperatures;

plt.plot(temperatures,m_var)
plt.show()