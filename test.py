import unittest
from mes import integrate_fun, integrate_fun_prime, B, l
from math import isclose

EPSILON = 1E-4

print("TESTS")

class TestStringMethods(unittest.TestCase):

  def test_E1E1_primes(self):
    result = 6
    I_result = integrate_fun_prime(1,1)[0]
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_E1E1(self):
    result = 2/9
    I_result = integrate_fun(1,1)[0]
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_E1E2_primes(self):
    result = -3
    I_result = integrate_fun_prime(1,2)[0]
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_E1E2(self):
    result = 1/18
    I_result = integrate_fun(1,2)[0]
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_B11(self):
    result = -6 -2/9
    I_result = B(1,1)
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_B12(self):
    result = 3 - 1/18
    I_result = B(1,2)
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_B21(self):
    result = 3 - 1/18
    I_result = B(2,1)
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_l1(self):
    result = 1/3
    I_result = l(1)
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

  def test_l2(self):
    result = 1/3
    I_result = l(2)
    self.assertIs(isclose(result, I_result, abs_tol=EPSILON), True)

if __name__ == '__main__':
    unittest.main()