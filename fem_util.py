from scipy.integrate import quad
from scipy.misc import derivative
import numpy as np
import numpy.linalg as la

def e_n(x, k, a=0, b=1, N=3):
  h=(b-a)/N
  x_k_prev = a + (k-1)*h
  x_k = a + k*h
  x_k_next = a + (k+1)*h
  if k == 0 and x >= a and x <= x_k_next:
    return (x_k_next - x)/h
  if k == N and x >= x_k_prev and x <= b:
    return (x - x_k_prev)/h
  elif x >= x_k_prev and x <= x_k:
    return (x - x_k_prev)/h
  elif x > x_k and x <= x_k_next:
    return (x_k_next - x)/h
  else:
    return 0

def B(f_u_v, u, v, a=0, b=1, N=3):
  return quad(f_u_v, a, b, args=(u,v))[0]

def l(f_x, v, a=0, b=1,):
  f_x_v = lambda x: f_x(x)*v(x)
  return quad(f_x_v, a, b, args=())[0]

def fem(f_x, f_u_v, u_shift=lambda x: 0, a=0, b=1, N=3, debug=False):
  d = 2 # matrix dimention
  matrix = np.zeros((d,d))
  e_base = [(lambda y: (lambda x: e_n(x, y+1)))(i) for i in range(d)]
  for n in range(d):
    for k in range (d):
      matrix[n][k] = B(f_u_v, e_base[k], e_base[n])

  vector = np.zeros((d,1))
  for i in range(d):
    vector[i] = l(f_x, e_base[i], a=a, b=b) - B(f_u_v, u_shift, e_base[i])

  if debug:
    print(matrix, "\n")
    print(vector, "\n")

  return la.solve(matrix, vector)
