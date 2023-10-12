#This file is to compute all the geodesics around Kerr spacetime.
#Jam Sadiq Aug 8th, 2016
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, cos, sin, pow, atan2
##================================================
#Initial Conditions
t0 = 0.0
x0 = 20.0   # geodesics starts along x-axis at this distance from BHs
y0 = 0.0
z0 = 0.0

tend = 2000000


rr = sqrt(x0*x0 + y0*y0 +z0*z0);
v0 =  1/sqrt(rr) #+ pert //kepler => pert =0 can change it -ve or +ve   
theta = atan2(y0,x0)
vx0 = -v0*sin(theta)
vy0 = v0*cos(theta)
vz0 = 0.0

##=================================================
#Importing metric from other file change name next line for import file

import gab_Binary as gab

get_gab = gab.get_gab
get_BHposition = gab.get_BHpositions
import RandGamma as geo
geo.set_gab(gab.get_gab)

geo.set_h(1.0e-3)

gab = get_gab(t0, x0, y0, z0) # where should I need to start coords?? Harmonic coords, Johnson Mcdaniel Paper
# Here mygab contains all metric components
x1, y1, z1, x2, y2, z2 = get_BHposition(t0)

print gab

test = geo.get_Riemann(t0, x0, y0, z0)

def mr(Svec):
    
    t = Svec[0]   #t,x,y,z
    x = Svec[1]
    y = Svec[2]
    z = Svec[3]
    u = Svec[4:8]       #ut,ux,uy,uz
    
    Gamma = geo.get_christoffel(t,x,y,z)
    out = np.ndarray(shape = 8, dtype = np.float64)
    out[0] = u[0]   #tdot
    out[1] = u[1]   #xdot
    out[2] = u[2]   #ydot
    out[3] = u[3]   #zdot
    
    for c in range(4):
        out[4+c] = 0  #setting zero and then filling in the geodesic equation
        
        for a in range(4):
            for b in range(4):
                out[4+c] -= Gamma[c, a, b]*u[a]*u[b]
    
    return out

def RK4(dt, yold, RHS):
    k1 = RHS(yold)
    ynew = yold + 0.5 * dt *k1
    k2 = RHS(ynew)
    ynew = yold + 0.5 * dt * k2
    k3 = RHS(ynew)
    ynew = yold + dt * k3
    k4 = RHS(ynew)
    ynew = yold + dt / 6 *(k1 + 2*k2 + 2*k3 + k4)
    return ynew


##==================================================
# Dot product
def gab_dot(v1, met, v2):
    v1down = np.dot(v1, met)
    return np.dot(v1down, v2)

#Uup = np.array((1,0,0,0), dtype=np.float64)
#vIphi = 1.0/(x0*sqrt(x0 -3))  #1./sqrt(x0)
Uup = np.array((0,vx0,vy0,vz0), dtype=np.float64) #change ydot =0.05 here
#making timelike four velocity by solving quadratic
a = gab[0,0]
b = 0
c = 0
for i in range(1,4):
  b = b + 2 *gab[0,i] * Uup[i]
  for j in range(1,4):
    c = c + gab[i,j]* Uup[i]*Uup[j]

c = c+1
d =np.sqrt(b**2 - 4*a*c)
Uup[0] = (-b - d)/(2*a)
print "ut0 ="
print Uup[0]
# is it correct??

#Udn = np.dot(gab, Uup)
#norm = np.dot(Uup, Udn)
#Uup = Uup / sqrt(-norm)

##=======================================================



yold  = np.ndarray(shape = 8, dtype = np.float64)   # i change shape=8 to 20
yold[0] = t0        #t0
yold[1] = x0     #x0
yold[2] = y0        #y0
yold[3] = z0        #z0
yold[4] = Uup[0]  # we need condition on it gab ua ub == -1
yold[5] = Uup[1]       # ux0
yold[6] = Uup[2]   #uy0  using gm/r = v**2/r, pow
yold[7] = Uup[3]    #uz0


def print_geodesic(outfile, time, y):
    R = geo.get_Riemann(y[0],y[1],y[2],y[3])
    u = y[4:8]
    Bx1, By1, Bz1, Bx2, By2, Bz2 = get_BHposition(y[0])
    xstar = y[1]#*np.cos(y[3])
    ystar = y[2]#*np.sin(y[3])
    zstar = y[3]
   #Bx1, By1, Bz1, Bx2, By2, Bz2 = get_BH_position(time) 

    
    outfile.write("{0:20.16e}".format(time))
    outfile.write(" {0:20.16e}".format(xstar))
    outfile.write(" {0:20.16e}".format(ystar))
    outfile.write(" {0:20.16e}".format(zstar))
    outfile.write(" {0:20.16e}".format(Bx1))
    outfile.write(" {0:20.16e}".format(By1))
    outfile.write(" {0:20.16e}".format(Bz1))
    outfile.write(" {0:20.16e}".format(Bx2))
    outfile.write(" {0:20.16e}".format(By2))
    outfile.write(" {0:20.16e}".format(Bz2))
    for i in range(8):
        outfile.write(" {0:20.16e}".format(y[i]))
    
    outfile.write("\n")


dt = 0.1  #step size which determine the accuracy of method increasing one decimal improve accracy by 4 time.
time = 0.0

#Defining automatic way to save the running output of the program

outfile=open("dataEinsteinBinay_x_{0}_v_{1}.txt".format(x0,v0), "w")

print_geodesic(outfile, time, yold)

x =  []
y =[]
xb1 = []
yb1 = []
xb2 = []
yb2 = []

for i in range(tend):
    ynew = RK4(dt, yold, mr)
    time = time + dt
    yold = ynew
    x1, y1, z1, x2, y2, z2 = get_BHposition(ynew[0])
    x.append(ynew[1])
    y.append(ynew[2])
    xb1.append(x1)
    yb1.append(y1)
    xb2.append(x2)
    yb2.append(y2)


    if (i%10 == 0):
        j = 1/10
        print_geodesic(outfile, time, ynew)
        #plt.plot(x,y)
        #plt.plot(xb1, yb1)
        #plt.plot(xb2, yb2)

        #plt.xlim(-x0+3, x0+3)
        #plt.xlim(-x0+3, x0+3)
        #plt.savefig('movie{0:08d}.png'.format(j))
        #plt.clf()
#The above if is to just output the required lines not the every iterations as those are unnecessary for analysis....


outfile.close()

plt.plot(x,y)
plt.plot(xb1, yb1) 
plt.plot(xb2, yb2)
plt.show()

