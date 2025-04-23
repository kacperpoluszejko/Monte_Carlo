import numpy as np
import matplotlib.pyplot as plt


# t, Dxx, Dyy, Dxy = np.loadtxt("Monte_carlo_6_10^2.txt", delimiter=" ", unpack=True, usecols=(0, 1, 2, 3))

# mask = np.isfinite(Dxx) & np.isfinite(Dyy) & np.isfinite(Dxy)
# Dxx = Dxx[mask]
# Dyy = Dyy[mask]
# Dxy = Dxy[mask]

# # Obliczenia
# N = len(Dxx)

# sumaxx = np.sum(Dxx)/N
# Dxx2 = Dxx*Dxx
# sumaxx2 = np.sum(Dxx2)/N

# sumayy = np.sum(Dyy)/N
# Dyy2 = Dyy*Dyy
# sumayy2 = np.sum(Dyy2)/N

# sumaxy = np.sum(Dxy)/N
# Dxy2 = Dxy*Dxy
# sumaxy2 = np.sum(Dxy2)/N

# print(sumaxx, sumayy, sumaxy)

# odchxx = np.sqrt((sumaxx2 - sumaxx*sumaxx)/N)
# odchyy = np.sqrt((sumayy2 - sumayy*sumayy)/N)
# odchxy = np.sqrt((sumaxy2 - sumaxy*sumaxy)/N)

# print(odchxx, odchyy, odchxy)

# N_table = [100, 1000, 10000, 100000]
# odchxx_table = [0.003075118454851255, 0.001527983906950655, 0.0012775425049057704, 0.0012646433933461711]
# odchyy_table = [0.0045630424964199666, 0.0017043555279236502, 0.0014093153689630023, 0.0012439086033146877]
# odchxy_table = [0.002465183782316827, 0.0009153619933331719, 0.0001418266634521953, 5.4040629589020794e-05]

# varxx_table = [0.8423184814814815, 0.9834586286286287, 1.0238962042042044, 1.005109978978979]
# varyy_table = [1.0216767377377378 , 0.9497157557557557, 1.0164626136136137, 1.013059781781782]
# varxy_table = [0.033665717888888884, -0.01452328849061061, 0.00817830850880881, -0.0036206448780780782]

# plt.errorbar(N_table, varxx_table, yerr=odchxx_table, fmt='o', capsize=5, markersize = 2, color = "black")
# plt.xlabel(r"$N_{max}$")
# plt.ylabel(r"$D_{xx}$")
# plt.xscale('log')
# plt.savefig("dxx.jpg")
# plt.show()

# plt.errorbar(N_table, varxy_table, yerr=odchxy_table, fmt='o', capsize=5, markersize = 2, color = "black")
# plt.xlabel(r"$N_{max}$")
# plt.ylabel(r"$D_{xy}$")
# plt.xscale('log')
# plt.savefig("dxy.jpg")
# plt.show()

# plt.errorbar(N_table, varyy_table, yerr=odchyy_table, fmt='o', capsize=5, markersize = 2, color = "black")
# plt.xlabel(r"$N_{max}$")
# plt.ylabel(r"$D_{yy}$")
# plt.xscale('log')
# plt.savefig("dyy.jpg")
# plt.show()

# plt.plot(t, Dxx, color="black", label = "Dxx")
# plt.plot(t, Dyy, color="red", label = "Dyy")
# plt.plot(t, Dxy, color="blue", label = "Dxy")
# plt.xlabel("t")
# plt.ylabel("D")
# plt.title(r"$N_{max} = 10^5$")
# plt.legend()
# plt.savefig("Monte_carlo_6_10^5")
# plt.show()

r = 0.5
R = 5  # promień
theta = np.linspace(0, 2 * np.pi, 100)
x = 3 + r * np.cos(theta)
y = r * np.sin(theta)
x1 =  R * np.cos(theta)
y1 =  R * np.sin(theta)

X1, Y1 = np.loadtxt("Monte_carlo_6_zad_2_t_1000.txt", delimiter=" ", unpack=True, usecols=(0, 1))
plt.figure(figsize=(6,6))
plt.plot(x,y,color = "blue")
plt.plot(x1,y1,color = "red")
plt.xlabel("x")
plt.ylabel("y")
plt.title("t=1000")
plt.scatter(X1, Y1, s=1, color='black')
plt.savefig("t=1000.jpg")
plt.show()

# X2, Y2 = np.loadtxt("Monte_carlo_6_zad_2_n_omega10_2.txt", delimiter=" ", unpack=True, usecols=(0, 1))
# plt.figure(figsize=(6,6))
# plt.xlabel("t")
# plt.ylabel("n")
# plt.title(r"$R_a = 0.1, \omega = 10$")
# plt.plot(X2, Y2, "red")
# plt.savefig("Ra_01_omega_10.png")
# plt.show()


# X1, Y1 = np.loadtxt("Monte_carlo_6_1.txt", delimiter=" ", unpack=True, usecols=(0, 1))
# X2, Y2 = np.loadtxt("Monte_carlo_6_2.txt", delimiter=" ", unpack=True, usecols=(0, 1))
# X3, Y3 = np.loadtxt("Monte_carlo_6_3.txt", delimiter=" ", unpack=True, usecols=(0, 1))

# plt.figure(figsize=(6,6))
# plt.scatter(X1, Y1, s=1, color='black', label = "t = 5")
# plt.scatter(X2, Y2, s=1, color='red', label = "t = 1")
# plt.scatter(X3, Y3, s=1, color='blue', label = "t = 0.1")
# plt.legend()
# plt.savefig("Monte_carlo_6_zad_1_1")
# plt.show()