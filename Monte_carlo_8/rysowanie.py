import numpy as np
import matplotlib.pyplot as plt
import math

X, Y = np.loadtxt("pcf_1.txt", delimiter=" ", unpack=True, usecols=(0, 1))

plt.plot(X,Y)
plt.show()