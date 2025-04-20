import matplotlib.pyplot as plt 
import numpy as np


data1 = np.loadtxt("monte_carlo_5_1.txt")
data2 = np.loadtxt("monte_carlo_5_2.txt")
data3 = np.loadtxt("monte_carlo_5_3.txt") #1000 losowań w systematycznej, 10^5 w warstwowej

counts1, bins1 = np.histogram(data1, bins=10)
counts2, bins2 = np.histogram(data2, bins=10)
counts3, bins3 = np.histogram(data3, bins=10)

bin_centers1 = 0.5 * (bins1[:-1] + bins1[1:])
bin_centers2 = 0.5 * (bins2[:-1] + bins2[1:])
bin_centers3 = 0.5 * (bins3[:-1] + bins3[1:])

X =np.linspace(-3, 3, 200)
Y =1 + np.tanh(X)

fig, ax1 = plt.subplots()


ax1.plot(bin_centers1, counts1, marker='o', linestyle='-', color='blue', label = "podstawowy")
ax1.plot(bin_centers2, counts2, marker='o', linestyle='-', color='black', label = "systematyczny")
ax1.plot(bin_centers3, counts3, marker='o', linestyle='-', color='red', label = "warstwowy")
ax1.set_xlabel("x")
ax1.set_ylabel(r"$N_m$")
plt.legend(loc = "lower center")

ax2 = ax1.twinx()
ax2.plot(X, Y,color = "green", label = r"g(x) = 1 + tanh(x)", linewidth = 0.5)
ax2.set_ylabel("g(x)")
plt.legend()
plt.title(r"$N = 10^3$, M=10")
plt.savefig("histogram3.png")
plt.show()


# # Dane
# x = [1, 2, 3, 4, 5]
# y1 = [10, 20, 25, 30, 40]  # Pierwsza skala (np. liczba sztuk)
# y2 = [100, 200, 300, 400, 500]  # Druga skala (np. złote)

# # Tworzenie pierwszego wykresu
# fig, ax1 = plt.subplots()

# color = 'tab:blue'
# ax1.set_xlabel('Dzień')
# ax1.set_ylabel('Liczba sztuk', color=color)
# ax1.plot(x, y1, color=color, label='Liczba sztuk')
# ax1.tick_params(axis='y', labelcolor=color)

# # Tworzenie drugiej osi Y
# ax2 = ax1.twinx()  # Użycie tej samej osi X
# color = 'tab:red'
# ax2.set_ylabel('Złote', color=color)
# ax2.plot(x, y2, color=color, linestyle='--', label='Kwota [zł]')
# ax2.tick_params(axis='y', labelcolor=color)

# # Dodanie tytułu i pokazanie wykresu
# plt.title('Podwójna skala Y')
# fig.tight_layout()  # Dopasowanie układu
# plt.show()