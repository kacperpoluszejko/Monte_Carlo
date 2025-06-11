import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("hist0.dat")
data2 = np.loadtxt("hist2.dat")

bin_centers = data[:, 0]           
f_numerical = data[:, 1]           
f_theoretical = data[:, 2]    

bin_centers2 = data2[:, 0]           
f_numerical2 = data2[:, 1]           
f_theoretical2 = data2[:, 2]       


plt.figure(figsize=(8, 5))
plt.plot(bin_centers, f_numerical, label="start", linewidth=2)
# plt.plot(bin_centers, f_theoretical, label="Rozkład teoretyczny", linestyle="--", linewidth=2)

plt.plot(bin_centers2, f_numerical2, label="end", linewidth=2)
# plt.plot(bin_centers2, f_theoretical2, label="Rozkład teoretyczny", linestyle="--", linewidth=2)

plt.xlabel("Prędkość [m/s]")
plt.ylabel("Rozkład cząstek")
# plt.title("it = 1000")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/zad_4_histogram.png")
plt.show()
