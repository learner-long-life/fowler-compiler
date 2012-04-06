from math import fabs, sqrt, sin, cos
from cmath import pi

class Matrix(object):
  def __init__(self, z11, z12, z21, z22):
    self.z11 = z11
    self.z12 = z12
    self.z21 = z21
    self.z22 = z22

def ca(a, b):
  return a + b

def cm(a, b):
  return a * b

def cc(a):
  return complex(a.real, -a.imag)

def md(M1, M2):
  z11 = ca(cm(cc(M1.z11), M2.z11), cm(cc(M1.z21), M2.z21));
  z22 = ca(cm(cc(M1.z12), M2.z12), cm(cc(M1.z22), M2.z22));

  return sqrt((2-abs(z11+z22))/2)

def pmd(a, b):
  return sqrt((2 - abs(sum(x*y for x, y in zip(a, b))))/2)

sqrt2o2 = sqrt(2)/2

G = Matrix(1, 0, 0, complex(cos(pi/6), sin(pi/6)))

I = Matrix(1, 0, 0, 1)
H = Matrix(sqrt2o2, sqrt2o2, sqrt2o2, -sqrt2o2)
X = Matrix(0, 1, 1, 0)
Z = Matrix(1, 0, 0, -1)

print md(I, G)
print md(H, G)
print md(X, G)
print md(Z, G)

print "----------"

PI = (1, 0, 0, 0)
PG = (cos(pi/12), 0, 0, -sin(pi/12))
print pmd(PI, PG)
