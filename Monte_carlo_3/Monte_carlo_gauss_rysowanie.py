import numpy as np
import matplotlib.pyplot as plt

# X, Y = np.loadtxt("Monte_carlo_gauss.txt", delimiter=" ", unpack=True, usecols=(0, 1))
# plt.figure(figsize=(6,6))
# plt.scatter(X, Y, s=1, color='black')
# plt.xlim(-5, 5)
# plt.ylim(-5, 5)
# plt.savefig("Carlo_1")
# plt.show()

# X2, Y2 = np.loadtxt("Monte_carlo_gauss_2.txt", delimiter=" ", unpack=True, usecols=(0, 1))
# print(X2)
# plt.figure(figsize=(6,6))
# plt.scatter(X2, Y2,  s=1, color='black')
# plt.xlim(-1.25, 1.25)
# plt.ylim(-1.25, 1.25)
# plt.savefig("Carlo_2")
# plt.show()

# X3, Y3 = np.loadtxt("Monte_carlo_gauss_3.txt", delimiter=" ", unpack=True, usecols=(0, 1))

# plt.figure(figsize=(6,6))
# plt.scatter(X3, Y3,  s=1, color='black')
# plt.xlim(-1.25, 1.25)
# plt.ylim(-1.25, 1.25)
# plt.savefig("Carlo_3")
# plt.show()

X4, Y4 = np.loadtxt("Monte_carlo_gauss_4.txt", delimiter=" ", unpack=True, usecols=(0, 1))
plt.figure(figsize=(6,6))
plt.scatter(X4, Y4,  s=1, color='black')
plt.xlim(-1.25, 1.25)
plt.ylim(-1.25, 1.25)
plt.savefig("Carlo_4")
plt.show()