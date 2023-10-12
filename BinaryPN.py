#This file is to compute all the geodesics around Kerr spacetime.
#Jam Sadiq Aug 8th, 2016
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, cos, sin, pow, atan2
import RK
plt.style.use("./presentation.mplstyle")
from mpl_toolkits.mplot3d import Axes3D



M = 1.0  #mass of BH
dt = 0.1 #for integration step
tend = 2000  #how long we want this run

#Initial Position and Velocity of orbiting body

#Initial Conditions  in  Cartesian x0 , y0, z0  

t0    = 0.0
x0 = 40.0   # geodesics starts along x-axis at this distance from the center of two BHs
y0 = 0.0
z0 = 0.0
#For precession we may need r and phi
rr  = sqrt(x0*x0 + y0*y0 +z0*z0);
phi = atan2(y0,x0)

def rphi(xv,yv,zv):
    Rr  = sqrt(xv*xv + yv*yv +zv*zv)
    Pphi = atan2(yv,xv)
    return Rr, Pphi


#initial velocities
pert = 0.0#10
v0 =  1.0/sqrt(rr)+ pert #kepler => pert =0 can change it -ve or +ve   
vx0 = v0*sin(phi)
vy0 = v0*cos(phi)
vz0 = 0.0






#Metric
import RandGamma as geo

from gab_Binary import get_gab as metric
from gab_Binary import get_BHpositions as Binary_position 

geo.set_gab(metric)                           #geo.set_gab(Kerr.get_gab)  #old 
#geo.set_gab(get_gabKerr)         
geo.set_h(1.0e-3)

#gab = metric(t0, R0, Thta0, Phi0) 
gab = metric(t0, x0, y0, z0) 


#For timelike geodesics
Uup = np.array((1,vx0, vy0, vz0), dtype=np.float64)
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
    Xx    = Svec[1]
    Yy = Svec[2]
    Zz = Svec[3]
    u = Svec[4:8]       #ut,ux,uy,uz
    
    Gamma = geo.get_christoffel(Tt,Xx,Yy,Zz)
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
yold[1] = x0        #x0
yold[2] = y0        #y0
yold[3] = z0        #z0
yold[4] = Uup[0]    # we need condition on it gab ua ub == -1
yold[5] = Uup[1]    #ux0
yold[6] = Uup[2]    #uy0  using gm/r = v**2/r, pow
yold[7] = Uup[3]    #uz0

#Data Saving for plots
BH1x = [Binary_position(t0)[0]]
BH1y = [Binary_position(t0)[1]]
BH1z = [Binary_position(t0)[2]]
BH2x = [Binary_position(t0)[3]]
BH2y = [Binary_position(t0)[4]]
BH2z = [Binary_position(t0)[5]]
TArray =[t0]
XArray=[x0]
YArray=[y0]
ZArray=[z0]

#print("t, x, y,z,  BH1x, BH1y,BH1z, BH2x, BH2y, BH2z")
#print(t, yold[1], yold[2], yold[3], Binary_position(t)[0], Binary_position(t)[1],Binary_position(t)[2], Binary_position(t)[3], Binary_position(t)[4], Binary_position(t)[5])

t = 0.0
while  t < tend:
  #(t, yold, dt) = RK.RK4_Step(t, yold, dt, RHS)#, tol=1.0e-10)
  (t, yold, dt) = RK.RK45_Step(t, yold, dt, RHS, tol=1.0e-12)
  #print(t, yold[1], yold[2], yold[3], Binary_position(t)[0], Binary_position(t)[1],Binary_position(t)[2], Binary_position(t)[3], Binary_position(t)[4], Binary_position(t)[5])
  BH1x.append(Binary_position(t)[0])
  BH1y.append(Binary_position(t)[1]) 
  BH1z.append(Binary_position(t)[2])
  BH2x.append(Binary_position(t)[3])
  BH2y.append(Binary_position(t)[4])
  BH2z.append(Binary_position(t)[5])
  TArray.append(t)
  XArray.append(yold[1])
  YArray.append(yold[2])
  ZArray.append(yold[3])
#  if (t%100==0.0):
#    print(XArray)
#    plt.plot(BH1x,BH1y, 'bo')
#    plt.plot(BH2x, BH2y, 'g*')
#    plt.plot(XArray, YArray, 'r--', label = "R ={0}".format(x0))
#    plt.legend()
#    plt.xlabel("X")
#    plt.ylabel("Y")
#    plt.show() 
#
plt.figure(figsize =(6,8))
plt.plot(BH1x,BH1y, 'bo')
plt.plot(BH2x, BH2y, 'g*')
plt.plot(XArray, YArray, 'r--')
plt.show()
#ThreeDplots
fig = plt.figure(figsize=(6,8))
ax = fig.gca(projection='3d')

ax.plot(BH1x[0],BH1y[0], BH1z[0],'k', marker='o', s = 12)
ax.plot(BH2x[0], BH2y[0], BH2z[0], 'k', marker= 'o', s = 12)

ax.plot(BH1x,BH1y, BH1z,'orange', marker='--')
ax.plot(BH2x, BH2y, BH2z, 'orange', marker= '--')
ax.plot(XArray, YArray,ZArray,'r')
ax.set_zlim(-0.1, 0.1)
plt.legend()
plt.title("Einsteinian case perturbed")
plt.show()

#save rphi data here for precession
