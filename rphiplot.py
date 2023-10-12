import numpy as np
import matplotlib.pyplot as plt
from math import *
import sys
data = np.loadtxt(sys.argv[1])
data = data.T
x = data[1]
y = data[2]
#z = data[3]
#r = np.sqrt(x**2 + y**2 + z**2)
r = np.sqrt(x**2 + y**2)
phi = np.unwrap(np.arctan2(y, x))
plt.plot(phi, r, linewidth=2)
plt.xlabel("phi")
plt.ylabel("r")
plt.title("Precession")
plt.show()

