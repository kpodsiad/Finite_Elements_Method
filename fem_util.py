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

def B(f_u_v, u, v, c_u_v=lambda u,v: 0, a=0, b=1, N=3):
  return quad(f_u_v, a, b, args=(u,v))[0] + c_u_v(u,v)

def l(f_x, v, f_v, a=0, b=1):
  f_x_v = lambda x: f_x(x)*v(x)
  return quad(f_x_v, a, b, args=())[0] + f_v(v)

def fem(f_x, f_u_v, u_shift=lambda x: 0, f_v=lambda v: 0, c_u_v=lambda u,v: 0, a=0, b=1, N=3, robin_left=False, robin_right=False, debug=False):
  d = N - 1 # matrix dimention
  if robin_left:
    d += 1
  if robin_right:
    d += 1

  # without left robin condition base looks like [e1,e2,e3...]
  # but with it is [e0,e1,e2,e3...]
  # so simply compute proper shift
  e_base_index_shift = 0 if robin_left else 1
  # must use another lambda because single lamda refer to the i via the reference
  e_base = [(lambda y: (lambda x: e_n(x, y+e_base_index_shift)))(i) for i in range(d)]

  matrix = np.array( [ [B(f_u_v, u=e_base[k], v=e_base[n], c_u_v=c_u_v) for k in range(d) ] for n in range(d) ] )

  vector = np.array( [ l(f_x=f_x, v=e_base[i], f_v=f_v, a=a, b=b) - B(f_u_v, u_shift, e_base[i])  for i in range(d) ] )

  if debug:
    print("matrix: \n", matrix, "\n")
    print("vector: \n",vector, "\n")

  return la.solve(matrix, vector)
