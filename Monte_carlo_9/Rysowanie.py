import numpy as np
import matplotlib.pyplot as plt

V = np.loadtxt('Monte_carlo_9_5.txt')  # domyślny separator to spacja
V2 = np.loadtxt('Monte_carlo_9_2.txt')
sigma = np.loadtxt('Monte_carlo_9_5_sigma.txt')
S = np.loadtxt('Monte_carlo_9_5_chains.txt')

plt.figure(figsize=(5, 5))
plt.imshow(V, origin='lower', cmap='seismic', vmin=-1, vmax=1, extent=[0, 20, 0, 20])
plt.colorbar(label='Potencjał')
plt.xlabel('i')
plt.ylabel('j')
plt.tight_layout()
plt.title(r"$V_{MC}$")
plt.savefig("V_MC_5.png")
plt.show()

plt.figure(figsize=(5, 5))
plt.imshow(np.abs(V-V2), origin='lower', cmap='seismic', vmin=-1, vmax=1, extent=[0, 20, 0, 20])
plt.colorbar()
plt.xlabel('i')
plt.ylabel('j')
plt.tight_layout()
plt.title(r"$|V_{MC} - V_{rel}|$")
plt.savefig("V_MCrel_5.png")
plt.show()

plt.figure(figsize=(5, 5))
plt.imshow(sigma, origin='lower', cmap='seismic', vmin=-1, vmax=1, extent=[0, 20, 0, 20])
plt.colorbar()
plt.xlabel('i')
plt.ylabel('j')
plt.tight_layout()
plt.title(r"$\sigma_V$")
plt.savefig("V_sigma_5.png")
plt.show()


plt.figure(figsize=(5, 5))
plt.imshow(S, origin='lower', cmap='seismic', vmin=-1, vmax=1, extent=[0, 20, 0, 20])
plt.colorbar()
plt.xlabel('i')
plt.ylabel('j')
plt.tight_layout()
plt.title(r"Absorbed chains")
plt.savefig("V_chains_5.png")
plt.show()