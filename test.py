import matplotlib.pyplot as plt
import json

gamma = [0.01, 1, 2, 3, 4]
f = 0.5
v = [(1 - f**(1/g))**g for g in gamma]
plt.plot(gamma, v)
plt.show()
