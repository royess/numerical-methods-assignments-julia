import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def runge(x):
    return 1 / (1 + x ** 2)


x_arr_plot = np.linspace(-5, 5, 500)
y_arr_plot = runge(x_arr_plot)

plt.figure(figsize=(20, 15))
plt.grid()

plt.plot(x_arr_plot, y_arr_plot, label='runge')

x1_arr = np.array([-5 + i for i in range(0,10+1)])
newton_poly = interpolate.lagrange(x1_arr, runge(x1_arr))
plt.plot(x_arr_plot, newton_poly(x_arr_plot), label='newton')

x2_arr = np.array([5 * np.cos((2 * i + 1) / 42 * np.pi) for i in range(0, 20+1)])
lagrange_poly = interpolate.lagrange(x2_arr, runge(x2_arr))
plt.plot(x_arr_plot, lagrange_poly(x_arr_plot), label='lagrange')

linear_interp = interpolate.interp1d(x1_arr, runge(x1_arr), 'linear')
plt.plot(x_arr_plot, linear_interp(x_arr_plot), label='linear')

hermite_interp = interpolate.PchipInterpolator(x1_arr, runge(x1_arr))
plt.plot(x_arr_plot, hermite_interp(x_arr_plot), label='hermite')

plt.legend()

plt.savefig('interpolations.jpg')