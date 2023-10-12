#This file is to compute all the geodesics around Kerr spacetime.
#Jam Sadiq Aug 8th, 2016
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, cos, sin, pow
import RK
plt.style.use("./presentation.mplstyle")
M = 1.0  #mass of BH
a = 0.0 #a = 0.0 => Schwarzschild ; for Kerr a can go 0 < a < 1
dt = 0.5 #for integration step
tend = 2000  #how long we want this run

#Initial Position and Velocity of orbiting body

#Initial Conditions  need same as Cartesian x0 = r , y0 = theta z0 = phi 
#We are using spherical coords here so we need translation to be in Cartesian

t0 = 0.0
R0 = 6.0        #geodesics at this distance from BHs 
Thta0 = np.pi/2.0  #mean we are on equator
Phi0 = 0.0        #We want velocity to body in phi direction

Rdot0 = 0.0 # vr
Thtadot0 = 0.0 #vtheta
Phidot0 = 1.0/(R0*sqrt(R0 -3))  # Worked out for circular orbit


#Metric
def get_Kerrgab(t, r, theta, phi):
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


#import gabKerr as Kerr
#get_gab = Kerr.get_gab  #old way if we get metric from gabKerr file  
import RandGamma as geo

geo.set_gab(get_Kerrgab)         
geo.set_h(1.0e-3)

gab = get_Kerrgab(t0, R0, Thta0, Phi0) 


#For timelike geodesics
Uup = np.array((1,Rdot0, Thtadot0, Phidot0), dtype=np.float64)
#For timelike geodesics
Aa = gab[0,0]
Bb = 0
Cc = 0
for i in range(1,4):
  Bb = Bb + 2 *gab[0,i] * Uup[i]
  for j in range(1,4):
    Cc = Cc + gab[i,j]* Uup[i]*Uup[j]

Cc = Cc+1
Dd =np.sqrt(Bb**2 - 4*Aa*Cc)
Uup[0] = (-Bb - Dd)/(2*Aa)
#print "check it is correct root +ve one for ut0 = "
print(Uup[0])





def RHS(time, Svec):
    Tt    = Svec[0]   #t,x,y,z
    Rr    = Svec[1]
    Tthta = Svec[2]
    Pphi  = Svec[3]
    u = Svec[4:8]       #ut,ux,uy,uz
    
    Gamma = geo.get_christoffel(Tt,Rr,Tthta,Pphi)
    out = np.ndarray(shape = 8, dtype = np.float64)
    out[0] = u[0]   #Tdot
    out[1] = u[1]   #Rdot
    out[2] = u[2]   #Thtadot
    out[3] = u[3]   #Phidot
    
    for c in range(4):
        out[4+c] = 0  #setting zero and then filling in the geodesic equation
        for a in range(4):
            for b in range(4):
                out[4+c] -= Gamma[c, a, b]*u[a]*u[b]
    return out

##======================================================
yold  = np.ndarray(shape = 8, dtype = np.float64)   # i change shape=8 to 20
yold[0] = t0        #t0
yold[1] = R0        #x0
yold[2] = Thta0        #y0
yold[3] = Phi0        #z0
yold[4] = Uup[0]    # we need condition on it gab ua ub == -1
yold[5] = Uup[1]    #ux0
yold[6] = Uup[2]    #uy0  using gm/r = v**2/r, pow
yold[7] = Uup[3]    #uz0

t = 0.0
TArray =[t]
XArray=[R0]
YArray=[Phi0]



while  t < tend:
  (t, yold, dt) = RK.RK4_Step(t, yold, dt, RHS)#, tol=1.0e-10)
  #(t, yold, dt) = RK.RK45_Step(t, yold, dt, RHS, tol=1.0e-10)
  TArray.append(t)
  XArray.append(yold[1])
  YArray.append(yold[3])
  
print(XArray)
plt.figure(figsize=(8,8))
plt.plot(XArray*np.cos(YArray), XArray*np.sin(YArray), 'r', label = "R = 6M, a={0}".format(a))
plt.legend(loc = 1, fontsize = 16)
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
