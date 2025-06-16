import numpy as np
import matplotlib.pyplot as plt

# Wczytaj dane pozycji i prędkości
data = np.loadtxt("rv.dat")

x = data[:, 0]
y = data[:, 1]

plt.figure(figsize=(5, 5))
plt.scatter(x, y, s=0.0001, color='black')
plt.xlabel("x")
plt.ylabel("y")
plt.title("koniec")
plt.xlim(0, 2)
plt.ylim(0, 0.5)
plt.gca().set_aspect('equal')
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/it_end.png")
plt.show()