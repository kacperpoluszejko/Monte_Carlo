import numpy as np
import matplotlib.pyplot as plt

#Zad 2 P_max = 100
t, x3, var = np.loadtxt("Monte_carlo_7_zad2.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2))
plt.errorbar(t, x3, yerr=var, fmt='o', capsize=3,elinewidth = 0.7,markeredgewidth=0.5,  markersize = 2, color = "black")
plt.xlabel("t")
plt.ylabel(r"$x3$")
plt.savefig("P_max=100.png")
plt.show()

#Zad 2 P_max = 5
t, x3, var = np.loadtxt("Monte_carlo_7_zad2_2.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2))
plt.errorbar(t, x3, yerr=var, fmt='o', capsize=3,elinewidth = 0.7,markeredgewidth=0.5,  markersize = 2, color = "black")
plt.xlabel("t")
plt.ylabel(r"$x3$")
plt.savefig("P_max=5.png")
plt.show()



#Zad 1
t, x1, x2, x3 = np.loadtxt("Monte_carlo_7_1.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2, 3))
plt.plot(t, x1, color = "black", label = "x1")
plt.plot(t, x2, color = "red", label = "x2")
plt.plot(t, x3, color = "blue", label = "x3")

# t, x1, x2, x3 = np.loadtxt("Monte_carlo_7_2.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2, 3))
# plt.plot(t, x1, color = "black")
# plt.plot(t, x2, color = "red")
# plt.plot(t, x3, color = "blue")

# t, x1, x2, x3 = np.loadtxt("Monte_carlo_7_3.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2, 3))
# plt.plot(t, x1, color = "black")
# plt.plot(t, x2, color = "red")
# plt.plot(t, x3, color = "blue")

# t, x1, x2, x3 = np.loadtxt("Monte_carlo_7_4.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2, 3))
# plt.plot(t, x1, color = "black")
# plt.plot(t, x2, color = "red")
# plt.plot(t, x3, color = "blue")

# t, x1, x2, x3 = np.loadtxt("Monte_carlo_7_5.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2, 3))
# plt.plot(t, x1, color = "black")
# plt.plot(t, x2, color = "red")
# plt.plot(t, x3, color = "blue")

plt.xlabel("t")
plt.ylabel("x1, x2, x3")
plt.legend()
plt.savefig("monte_carlo_7_1.png")
plt.show()