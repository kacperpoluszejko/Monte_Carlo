import numpy as np
import matplotlib.pyplot as plt

# Wczytaj dane pozycji i prędkości
data = np.loadtxt("wynikirv_1.dat")

x = data[:, 0]
y = data[:, 1]

plt.figure(figsize=(5, 5))
plt.scatter(x, y, s=0.01, color='black')
plt.xlabel("x")
plt.ylabel("y")
plt.title("iteracja 1")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal')
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/it_1.png")
plt.show()