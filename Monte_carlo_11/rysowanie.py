import numpy as np
import matplotlib.pyplot as plt
import matplotlib

data = np.loadtxt('Monte_carlo_11.txt')  
plt.imshow(data, norm = matplotlib.colors.LogNorm())
plt.colorbar()
plt.show()