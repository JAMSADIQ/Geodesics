import numpy as np
""" 
A very simple RK4 and RK45 implementation and vecrtor and scalar tests. 
"""


default_tol = 1.0e-10

def RK4_Step(xold, yold, dx, fxy):
  k1 = fxy(xold, yold)
  k2 = fxy(xold+ 0.5 * dx, yold + 0.5 * dx * k1)
  k3 = fxy(xold+ 0.5 * dx, yold + 0.5 * dx * k2)
  k4 = fxy(xold+ 1.0 * dx, yold + 1.0 * dx * k3)

  xnew = xold + dx
  ynew = yold + (dx / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
  return xnew, ynew, dx



# Adaptive RK45
def RK45_Step(xold, yold, dx, fxy, tol=default_tol):
  k1 = dx * fxy(xold, yold)
  k2 = dx * fxy(xold + 1.0/4.0 * dx, yold + 1.0/4.0 * k1)
  k3 = dx * fxy(xold + 3.0/8.0 * dx, yold + 3.0/32.0 * k1 + 9.0/32.0 * k2)
  k4 = dx * fxy(xold + 12.0/13.0 * dx, yold + 1932.0/2197.0 * k1 -\
                 7200.0/2197.0 * k2 + 7296.0/2197.0 * k3)
  k5 = dx * fxy(xold + dx, yold + 439.0/216.0 * k1 - 8.0 * k2 +\
                 3680.0/513.0 * k3 - 845.0/4104.0 * k4)
  k6 = dx * fxy(xold + 1.0/2.0 * dx, yold -8.0/27.0 * k1 + 2.0 * k2 -\
                 3544.0/2565.0 * k3 + 1859.0/4104.0 * k4 - 11.0/40.0 * k5)

  # higher order accuracy ynew
  ynew = yold + 16.0/135.0 * k1 + 6656.0/12825.0 * k3 +\
             28561.0/56430.0 * k4 -9.0/50.0 * k5 + 2.0/55.0 * k6

  # lower order accuracy ynew
  ynew4 = yold + 25.0/216.0 * k1 + 1408.0/2565.0 * k3 +\
             2197.0/4104.0 * k4 -1.0/5.0 * k5

  xnew = xold + dx

  err= abs(np.array(ynew-ynew4)).max()

#  if err < tol * 1.0e-3:
#    dtc = 2* dt
#  else:

  dxnew = 0.9 * dx * (tol/err)**(0.25)

  if err > tol:
    ynew = yold
    xnew = xold

  return (xnew, ynew, dxnew)


def rhs(x, y):
  return x**6 - y


def exact(x):
  res = 720 - 720*np.exp(-x)
  res += x*(-720 + x*(360 + x*(-120 + x*(30 + (-6 + x)*x))))
  return res


def test(xend, dx):
  y = 0
  x = 0
  while  x < xend:
    (x, y, dx) = RK4_Step(x, y, dx, rhs)
    x = x + dx
    print(x, y - exact(x))

def test(xend, dx):
  y = 0
  x = 0
  while  x < xend:
    (x, y, dx) = RK45_Step(x, y, dx, rhs, tol=1.e-13)
    print(x, y - exact(x))


def rhs2(x,y):
  return np.array((-y[1], y[0]))

def exact2(x):
  return np.array((-np.sin(x), np.cos(x)))

def test2(xend, dx):
  y = np.array((0,1), dtype=np.float64)
  x = 0
  while  x < xend:
    (x, y, dx) = RK4_Step(x, y, dx, rhs2)
    print(x, np.linalg.norm(y-exact2(x)))

