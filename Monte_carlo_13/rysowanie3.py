import numpy as np
import matplotlib.pyplot as plt

# Wczytaj dane pozycji i prędkości
data1 = np.loadtxt("wyniki/nptv_1000.dat")
data2 = np.loadtxt("wyniki/nptv_1250.dat")
data3 = np.loadtxt("wyniki/nptv_1500.dat")
data4 = np.loadtxt("wyniki/nptv_2000.dat")

x1 = data1[:, 0]
n1 = data1[:, 1]
p1 = data1[:, 2]
T1 = data1[:, 3]
v1 = data1[:, 4]

x2 = data2[:, 0]
n2 = data2[:, 1]
p2 = data2[:, 2]
T2 = data2[:, 3]
v2 = data2[:, 4]

x3 = data3[:, 0]
n3 = data3[:, 1]
p3 = data3[:, 2]
T3 = data3[:, 3]
v3 = data3[:, 4]

x4 = data4[:, 0]
n4 = data4[:, 1]
p4 = data4[:, 2]
T4 = data4[:, 3]
v4 = data4[:, 4]

plt.plot(x1, n1, label = "it 1000")
plt.plot(x2, n2, label = "it 1250")
plt.plot(x3, n3, label = "it 1500")
plt.plot(x4, n4, label = "it 2000")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.title("gęstość")
plt.savefig("plots/n_3")
plt.show()

plt.plot(x1, p1, label = "it 1000")
plt.plot(x2, p2, label = "it 1250")
plt.plot(x3, p3, label = "it 1500")
plt.plot(x4, p4, label = "it 2000")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.title("ciśnienie")
plt.savefig("plots/p_3")
plt.show()

plt.plot(x1, T1, label = "it 1000")
plt.plot(x2, T2, label = "it 1250")
plt.plot(x3, T3, label = "it 1500")
plt.plot(x4, T4, label = "it 2000")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.title("temperatura")
plt.savefig("plots/T_3")
plt.show()

plt.plot(x1, v1, label = "it 1000")
plt.plot(x2, v2, label = "it 1250")
plt.plot(x3, v3, label = "it 1500")
plt.plot(x4, v4, label = "it 2000")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.title("prędkość")
plt.savefig("plots/v_3")
plt.show()
