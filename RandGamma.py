#This file is to use the finite differencing to compute derivatives of metric tesor to compute the christofell symbols and Riemann tensor components. This file will be imported in main file where we will compute and save all the geodesics and scalars.
#geometry.py
#JamSadiq August 8th, 2016

import numpy as np
h_finite_diff = 1.0e-3
#I need to change h_finite_diff to make improvements, but there is a limit how much h can be change.

#=====================================================
get_gab = None

def set_gab(gab):
    global get_gab
    get_gab = gab


def set_h(h):
    global h_finite_diff
    h_finite_diff = h

#=======================================================
#def tderiv(f,t,x,y,z):
#   h = h_finite_diff
#   return (f(t+h,x,y,z)-f(t-h,x,y,z))*0.5/h
#def xderiv(f,t,x,y,z):
#    h = h_finite_diff
#    return (f(t,x+h,y,z)-f(t,x-h,y,z))*0.5/h
#def yderiv(f,t,x,y,z):
#    h = h_finite_diff
#    return (f(t,x,y+h,z)-f(t,x,y-h,z))*0.5/h
#def zderiv(f,t,x,y,z):
#    h = h_finite_diff
#    return (f(t,x,y,z+h)-f(t,x,y,z-h))*0.5/h
#=========================================================


#Extra Higher order derivatives

#=======================================================
def tderiv(f,t,x,y,z):
    h = h_finite_diff
    return (-f(t+2*h,x,y,z)+ 8*f(t+h,x,y,z)-8*f(t-h,x,y,z)+f(t-2*h,x,y,z))/(12*h)
def xderiv(f,t,x,y,z):
    h = h_finite_diff
    return (-f(t,x+2*h,y,z)+ 8*f(t,x+h,y,z)-8*f(t,x-h,y,z)+f(t,x-2*h,y,z))/(12*h)
def yderiv(f,t,x,y,z):
    h = h_finite_diff
    return  (-f(t,x,y+2*h,z)+ 8*f(t,x,y+h,z)-8*f(t,x,y-h,z)+f(t,x,y-2*h,z))/(12*h)
def zderiv(f,t,x,y,z):
    h = h_finite_diff
    return  (-f(t,x,y,z+2*h)+ 8*f(t,x,y,z+h)-8*f(t,x,y,z-h)+f(t,x,y,z-2*h))/(12*h)
#=========================================================

def get_christoffel(t,x,y,z):
    t=np.float64(t)
    x =np.float64(x)
    y = np.float64(y)
    z = np.float64(z)
    
    gab = get_gab(t,x,y,z)
    inverse_gab =  np.linalg.inv(gab)
    
    
    allgabderiv=[0,0,0,0]
    allgabderiv[0] = tderiv(get_gab, t,x,y,z )
    allgabderiv[1] = xderiv(get_gab, t,x,y,z )
    allgabderiv[2] = yderiv(get_gab, t,x,y,z )
    allgabderiv[3] = zderiv(get_gab, t,x,y,z )
    
    christoffel = np.ndarray(shape=(4,4,4), dtype=np.float64)
    
    for c in range(4):
        for a in range(4):
            for b in range(4):
                christoffel[c, a, b] = 0
                for d in range(4):
                    christoffel[c,a,b] += 0.5 * inverse_gab[c, d]*\
                        (allgabderiv[a][b,d] + allgabderiv[b][a,d] -\
                         allgabderiv[d][a,b])
    return christoffel

#===================================================================
#Riemann Tensor Computation


def get_Riemann(t,x,y,z):
    t=np.float64(t)
    r =np.float64(x)
    theta = np.float64(y)
    phi = np.float64(z)
    
    gab = get_gab(t,x,y,z)
    inverse_gab =  np.linalg.inv(gab)
    Gamma = get_christoffel(t,x,y,z)
    
    
    allgammaderiv=[0,0,0,0]
    allgammaderiv[0] = tderiv(get_christoffel, t,x,y,z )
    allgammaderiv[1] = xderiv(get_christoffel, t,x,y,z )
    allgammaderiv[2] = yderiv(get_christoffel, t,x,y,z )
    allgammaderiv[3] = zderiv(get_christoffel, t,x,y,z )
    
    Riemann = np.ndarray(shape=(4,4,4,4), dtype=np.float64)
    
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    Riemann[a,b,c,d] = 0
                    for e in range(4):
                        Riemann[a,b,c,d] += gab[a, e]*(allgammaderiv[c][e,b,d] - allgammaderiv[d][e,b,c])
                        for k in range(4):
                            Riemann[a,b,c,d] +=gab[a,e]* (Gamma[e,c,k]*Gamma[k,b,d]-   Gamma[e,d,k]*Gamma[k,b,c])
    
    return Riemann




