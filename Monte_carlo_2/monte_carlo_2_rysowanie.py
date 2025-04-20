import matplotlib.pyplot as plt 
import numpy as np


data = np.loadtxt("monte_carlo_2_1.txt")
print(len(data))
X = np.linspace(0,1,10000)
y = 0.8*(1+X-X*X*X)
plt.plot(X,y, color = "black")
plt.hist(data, bins=10, edgecolor='black', color='royalblue', density = True)
plt.ylim(0.8, 1.15)
plt.title("Metoda eliminacji")
plt.savefig("Metoda_eliminacji.jpg")
plt.show()

def dystrybuanta(x):
     return 0.8*(x + ((x**2)/2) - ((x**4)/4))

def p_i(i):
     x1 = i*0.1
     x2 = (i+1)*0.1
     return dystrybuanta(x2) - dystrybuanta(x1)


liczby, przedzialy = np.histogram(data, bins=10)
N = sum(liczby)
chi = 0

for i in range (len(liczby)):
    chi = chi + (liczby[i] - p_i(i)*N)**2/(p_i(i)*N)

    print(i,"  ", p_i(i))

print(chi)