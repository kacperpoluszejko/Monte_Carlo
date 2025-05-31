import numpy as np
import matplotlib.pyplot as plt

u, x = np.loadtxt("Monte_carlo_10_u4_3.txt", delimiter=" ", unpack=True, usecols=(0, 1))  # domyślny separator to spacja
u_teo, x_teo = np.loadtxt("Monte_carlo_10_u4_teo.txt", delimiter=" ", unpack=True, usecols=(0, 1))  # domyślny separator to spacja

plt.scatter(x,u, color = "red", s = 5, label = "MC")
plt.plot(x_teo, u_teo, label = "exact")
plt.title(r"$n_{paths} = 100000, t = 35$")
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.savefig("t=35_n=100000.png")
plt.show()


plt.plot(x_teo, np.abs(u_teo - u), label = "exact")
plt.title(r"$n_{paths} = 100000, t = 35$")
plt.xlabel("x")
plt.ylabel(r"$|u-u_{teo}|$")
plt.ylim(0, 0.002)
plt.savefig("t=35_n=100000_eroor.png")
plt.show()
