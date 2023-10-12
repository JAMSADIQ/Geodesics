import numpy as np
from math import sqrt, pow, cos, sin
M = 1.0
a = 0.8 #0.5
def get_gab(t, r, theta, phi):
  r2 = r*r
  a2 = a*a
  rho2 = r2 + a2*cos(theta)*cos(theta)
  Del  = r2 -2*M*r +a2
  Sigma = (r2 +a2)* (r2 +a2)  - a2*Del*sin(theta)*sin(theta)
  g00 =  -1+ (2*M*r)/rho2 
  g01 = 0
  g02 = 0
  g03 = -2*M*a*r*sin(theta)*sin(theta)/rho2 
  g11 = rho2/Del
  g12 = 0
  g13 = 0 
  g22 = rho2
  g23 = 0
  g33 = Sigma*sin(theta)*sin(theta)/rho2

  return np.array(((g00, g01, g02, g03), (g01, g11, g12, g13), (g02, g12, g22, g23),(g03, g13, g23, g33)), dtype=np.float64)
