import RK
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


binary_separation  = 100.0

mass_ratio = 1.0
M = 1.0
G = 1.0


def Kepler(t, y):

  """ Implements  d\vec X/dt = \vec V, d\vec V/dt = - M  G  \vec X/|X|^3.
  Takes as input a time and
  a 6-d vector. The first three components are the components of \vec
  X, the next three are the components of \vec V. The parameter "t" is
  not used.
  """

  position = y[0:3]
  velocity = y[3:6]

  # position_dot = velocity
  pos_dot = velocity;

  r = np.sqrt(np.dot(position,position))
  r3 = r*r*r

  # formula for velocity_dot
  vel_dot = - position  * M * G / r3;

  # y_dot 
  return np.concatenate((pos_dot, vel_dot))

def Binary_position(t):

  """ Gives x,y,z coordinated of the two components of a binary in a circular orbit. 
  """

  binary_frequency = np.sqrt(M/binary_separation**3)
  binary_angle = t * binary_frequency
  X1 = binary_separation * np.array([np.cos(binary_angle), np.sin(binary_angle),0]) / (1.0 + mass_ratio)
  X2 = -X1 
  return X1, X2
  

def Binary(t, y):

  """ Source term for a binary. It has components m1 and m2. Currently
  assumes circular orbits. Can adjust binary separation and mass
  ratio.
  """

  m1 = M * mass_ratio / (1.+mass_ratio)
  m2 = M / (1.+mass_ratio)

  X1, X2 = Binary_position(t)

  position = y[0:3]
  velocity = y[3:6]

  # position_dot = velocity
  pos_dot = velocity;

  p1 = position - X1
  p2 = position - X2

  r1 = np.sqrt(np.dot(p1,p1))
  r2 = np.sqrt(np.dot(p2,p2))

  # formula for velocity_dot
  vel_dot = - p1  * m1 * G / r1**3 - p2 * m2 * G / r2**3

  # y_dot 
  return np.concatenate((pos_dot, vel_dot))



# Initial position
x0 = 400 
y0 = 0
z0 = 0

#Initial velocity
Vx =0
Vy =1*np.sqrt(1.0/x0)
Vz =0

tend = 200000
dt = 1.0
y = np.array((x0,y0,z0, Vx, Vy, Vz), dtype=np.float64)
t = 0

#TimeArray=[t]
XArray=[x0]
YArray=[y0]

def plot_output(t, y, i):
    Xx, Xy, Xz, Vx, Vy, Vz = y
    plt.plot(y[0], y[1],  color='green', marker='o', linestyle='dashed', markersize=5)
    BH1, BH2 = Binary_position(t)
    plt.plot(BH1[0], BH1[1], color='black' , marker='o', linestyle='dashed', markersize=4)
    plt.plot(BH2[0], BH2[1], color='orange' , marker='o', linestyle='dashed', markersize=4)


    plt.xlim(-415,415)
    plt.ylim(-415,415)

    #plt.savefig('movie{0:05d}.png'.format(i))  #commented 
    #added black hole positions in the output 
    print t, y[0], y[1], y[2], BH1[0], BH1[1],BH1[2], BH2[0], BH2[1], BH2[2], y[3], y[4], y[5]  
#commented

    plt.clf()


i=0

delta_t_dump =1
t_dump = 0

plot_output(t, y, i)


while  t < tend:
  told = t
 #(t, y, dt) = RK.RK45_Step(t, y, dt, Kepler, tol=1.0e-12)
  (t, y, dt) = RK.RK45_Step(t, y, dt, Binary, tol=1.0e-12)
 
  if (t>= t_dump + delta_t_dump):

    i = i + 1
    plot_output(t, y, i)
    t_dump = t




#plt.plot(XArray, YArray)

#plt.axes().set_aspect('equal','datalim')
#plt.show()

