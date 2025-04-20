import numpy as np
import matplotlib.pyplot as plt

# X1, Y1 = np.loadtxt("Monte_carlo_4_1.txt", delimiter=" ", unpack=True, usecols=(0, 1))
# X2, Y2 = np.loadtxt("Monte_carlo_4_2.txt", delimiter=" ", unpack=True, usecols=(0, 1))
# plt.figure(figsize=(6,6))
# plt.scatter(X1, Y1, s=1, color='black')
# plt.scatter(X2, Y2, s=1, color='blue')
# plt.xlim(-4, 8)
# plt.ylim(-6, 6)
# plt.savefig("Carlo_4_1")
# plt.show()

X, Y, Z = np.loadtxt("Monte_carlo_4_4.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2))
X1, Y1, Z1 = np.loadtxt("Monte_carlo_4_3.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2))

plt.errorbar(Z1, X1, yerr=Y1, fmt='o', color='red', label=r'$\alpha=A$', capsize=4)
plt.xscale('log')
plt.errorbar(Z, X, yerr=Y, fmt='o', color='black', label=r'$\alpha=B$', capsize=4)
plt.xscale('log')
plt.title(r"$x_a = 0$")
plt.legend()
plt.savefig("Carlo_4_2.jpg")
plt.show()