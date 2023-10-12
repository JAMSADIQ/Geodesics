import numpy as np
from math import *
import sys

in_filename=sys.argv[1]
out_name=sys.argv[1]+"r_versus_phi"
data = np.loadtxt(in_filename)
data = data.T
x = data[1]
y = data[2]
#z = data[3]
#r = np.sqrt(x**2 + y**2 + z**2)
r = np.sqrt(x**2 + y**2)
phi = np.unwrap(np.arctan2(y, x))
np.savetxt(out_name, np.stack((phi,r), axis=-1))

