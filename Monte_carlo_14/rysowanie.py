import numpy as np
import matplotlib.pyplot as plt


energy = np.loadtxt("energia.txt")
variance = np.loadtxt("wariancja.txt")



plt.imshow(variance, origin = "lower")
plt.colorbar()
plt.show()

variance = np.log(variance + 10**(-20))

plt.imshow(variance, origin = "lower")
plt.colorbar()
plt.show()