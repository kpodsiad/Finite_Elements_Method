from scipy.integrate import quad
from scipy.misc import derivative
import numpy as np
import numpy.linalg as la

# given equation u'' - u' = 1, domain=[0,1], u(0)=0, u(1)=0,

a=0
b=1
n=3
h=(b-a)/n

def e_n(x, k):
  x_k_prev = a + (k-1)*h
  x_k = a + k*h
  x_k_next = a + (k+1)*h
  if x >= x_k_prev and x <= x_k:
    return (x - x_k_prev)/h
  elif x > x_k and x <= x_k_next:
    return (x_k_next - x)/h
  else:
    return 0

def fun(x, n=1, k=1):
  return e_n(x, n)*e_n(x, k)

def fun_prime(x, n=1, k=1):
  return derivative(e_n, x, dx=1e-6, args=(n,)) * derivative(e_n, x, dx=1e-6, args=(k,))

def integrate_fun(n=1, k=1):
  return quad(fun, a, b, args=(n,k))

def integrate_fun_prime(n=1, k=1):
  return quad(fun_prime, a, b, args=(n,k))

def B(n=1, k=1):
  return -integrate_fun_prime(n,k)[0] - integrate_fun(n,k)[0]

def l(n=1):
  return quad(e_n, a, b, args=(n,))[0]

def fem():
  d = 2 # matrix dimention
  matrix = np.zeros((d,d))
  for k in range(d):
    for n in range (d):
      val = B(n+1, k+1)
      matrix[n][k] = val

  vector = np.zeros((d,1))
  for i in range(d):
    vector[i] = l(i+1)

  return la.solve(matrix, vector)
