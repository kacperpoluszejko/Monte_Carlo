import numpy as np
import matplotlib.pyplot as plt

# Wczytaj dane pozycji i prędkości
data = np.loadtxt("wyniki/nptv_10.dat")

x = data[:, 0]
p = data[:, 2]
T = data[:, 3]

# Główna oś (lewa oś y)
fig, ax1 = plt.subplots()

ax1.plot(x, T, 'b-', label='T')
ax1.set_xlabel('x')
ax1.set_ylabel('T', color='r')
ax1.tick_params(axis='y', labelcolor='r')

# Druga oś y (po prawej stronie)
ax2 = ax1.twinx()
ax2.plot(x, p, 'r--', label='p')
ax2.set_ylabel('p', color='b')
ax2.tick_params(axis='y', labelcolor='b')

# Tytuł i siatka
plt.title('it 10')
ax1.grid(True)

# Wyświetlenie
plt.savefig("plots/zad_4_it_10")
plt.show()