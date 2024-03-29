import unittest
from scipy.misc import derivative
from fem_util import B, e_n, l, fem
from math import isclose

EPSILON = 1E-6
a=0
b=1
e0 = lambda x: e_n(x, 0)
e1 = lambda x: e_n(x, 1)
e2 = lambda x: e_n(x, 2)
e3 = lambda x: e_n(x, 3)

case1_f_x = lambda x: 1
case1_f_u_v = lambda x,u,v: -derivative(u, x, dx=1e-6, args=())*derivative(v, x, dx=1e-6, args=()) - u(x)*v(x)

case2_f_x = lambda x: x
case2_f_u_v = lambda x,u,v: -derivative(u, x, dx=1e-6, args=())*derivative(v, x, dx=1e-6, args=()) + derivative(u, x, dx=1e-6, args=())*v(x)

case3_f_x = lambda x: x
case3_f_u_v = lambda x,u,v: derivative(u, x, dx=1e-6, args=())*derivative(v, x, dx=1e-6, args=())
case3_f_v = lambda v: -v(a)

case4_f_x = lambda x: 1
case4_f_u_v = lambda x,u,v: -derivative(u, x, dx=1e-6, args=())*derivative(v, x, dx=1e-6, args=()) + u(x)*v(x)
case4_f_v = lambda v: +v(b)
case4_c_u_v = lambda u,v: u(b)*v(b)
case4_u_shift = lambda x: 2*e0(x)


print("TESTS")

class TestFemMethods(unittest.TestCase):

  def test_case1_E1E1_primes(self):
    f = lambda x,u,v: -derivative(u, x, dx=1e-6, args=())*derivative(v, x, dx=1e-6, args=())
    result = isclose(B(f, u=e1, v=e1), -6, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_E1E1(self):
    f = lambda x,u,v: - u(x)*v(x)
    result = isclose(B(f, u=e1, v=e1), -2/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_E1E2_primes(self):
    f = lambda x,u,v: -derivative(u, x, dx=1e-6, args=())*derivative(v, x, dx=1e-6, args=())
    result = isclose(B(f, u=e1, v=e2), 3, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_E1E2(self):
    f = lambda x,u,v: - u(x)*v(x)
    result = isclose(B(f, u=e1, v=e2), -1/18, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_B11(self):
    result = isclose(B(case1_f_u_v, u=e1, v=e1), -6 -2/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_B12(self):
    result = isclose(B(case1_f_u_v, u=e1, v=e2), 3 - 1/18, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_B21(self):
    result = isclose(B(case1_f_u_v, u=e2, v=e1), 3 - 1/18, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_B22(self):
    result = isclose(B(case1_f_u_v, u=e2, v=e2), -6 -2/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_l1(self):
    result = isclose(l(case1_f_x, e1, f_v=lambda v: 0), 1/3, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_l2(self):
    result = isclose(l(case1_f_x, e2, f_v=lambda v: 0), 1/3, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case1_fem(self):
    result = fem(case1_f_x, case1_f_u_v)
    self.assertIs(isclose(result[0], -6/59, abs_tol=EPSILON), True)
    self.assertIs(isclose(result[1], -6/59, abs_tol=EPSILON), True)

  def test_case2_B11(self):
    result = isclose(B(case2_f_u_v, u=e1, v=e1), -6, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case2_B12(self):
    result = isclose(B(case2_f_u_v, u=e1, v=e2), 2.5, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case2_B21(self):
    result = isclose(B(case2_f_u_v, u=e2, v=e1), 3.5, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case2_B22(self):
    result = isclose(B(case2_f_u_v, u=e2, v=e2), -6, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case2_l1(self):
    val = l(case2_f_x, e1, f_v=lambda v: 0) - B(case2_f_u_v, u=e3, v=e1)
    result = isclose(val, 1/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case2_l2(self):
    val = l(case2_f_x, e2, f_v=lambda v: 0) - B(case2_f_u_v, u=e3, v=e2)
    result = isclose(val, 2/9 - 3.5, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case2_fem(self):
    result = fem(case2_f_x, case2_f_u_v, u_shift=e3)
    self.assertIs(isclose(result[0], 0.396534, abs_tol=EPSILON), True)
    self.assertIs(isclose(result[1], 0.711519, abs_tol=EPSILON), True)

  def test_case3_B00(self):
    result = isclose(B(case3_f_u_v, u=e0, v=e0), 3, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case2_B10(self):
    val = B(case3_f_u_v, u=e1, v=e0)
    result = isclose(val, -3, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_B20(self):
    result = isclose(B(case3_f_u_v, u=e2, v=e0), 0, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_B11(self):
    result = isclose(B(case3_f_u_v, u=e1, v=e1), 6, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_B12(self):
    result = isclose(B(case3_f_u_v, u=e1, v=e2), -3, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_B22(self):
    result = isclose(B(case3_f_u_v, u=e2, v=e2), 6, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_l0(self):
    val = l(case3_f_x, e0, f_v=case3_f_v)
    result = isclose(val, 1/54 - 1, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_l1(self):
    val = l(case3_f_x, e1, f_v=case3_f_v)
    result = isclose(val, 1/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_l2(self):
    val = l(case3_f_x, e2, f_v=case3_f_v)
    result = isclose(val, 2/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case3_fem(self):
    result = fem(f_x=case3_f_x, f_u_v=case3_f_u_v, f_v=case3_f_v, robin_left=True)
    self.assertIs(isclose(result[0], -5/6, abs_tol=EPSILON), True)
    self.assertIs(isclose(result[1], -41/81, abs_tol=EPSILON), True)
    self.assertIs(isclose(result[2], -35/162, abs_tol=EPSILON), True)

  def test_case4_B11(self):
    val = B(case4_f_u_v, u=e1, v=e1, c_u_v=case4_c_u_v)
    result = isclose(val, -6 + 2/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_B21(self):
    result = isclose(B(case4_f_u_v, u=e2, v=e1, c_u_v=case4_c_u_v), 3 + 1/18, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_B31(self):
    result = isclose(B(case4_f_u_v, u=e3, v=e1, c_u_v=case4_c_u_v), 0, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_B22(self):
    result = isclose(B(case4_f_u_v, u=e2, v=e2, c_u_v=case4_c_u_v), -6 + 2/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_B32(self):
    result = isclose(B(case4_f_u_v, u=e3, v=e2, c_u_v=case4_c_u_v), 3 + 1/18, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_B33(self):
    val = B(case4_f_u_v, u=e3, v=e3, c_u_v=case4_c_u_v)
    result = isclose(val, -2 + 1/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_l1(self):
    val = l(case4_f_x, e1, f_v=case4_f_v) - B(case4_f_u_v, u=case4_u_shift, v=e1)
    result = isclose(val, 1/3 - 6 - 1/9, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_l2(self):
    val = l(case4_f_x, e2, f_v=case4_f_v) - B(case4_f_u_v, u=case4_u_shift, v=e2)
    result = isclose(val, 1/3, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_l3(self):
    val = l(case4_f_x, e3, f_v=case4_f_v) - B(case4_f_u_v, u=case4_u_shift, v=e3)
    result = isclose(val, 1/6 + 1, abs_tol=EPSILON)
    self.assertIs(result, True)

  def test_case4_fem(self):
    result = fem(f_x=case4_f_x, f_u_v=case4_f_u_v, f_v=case4_f_v, c_u_v=case4_c_u_v, u_shift=case4_u_shift, robin_right=True)
    self.assertIs(isclose(result[0], 21601/49706, abs_tol=EPSILON), True)
    self.assertIs(isclose(result[1], -26572/24853, abs_tol=EPSILON), True)
    self.assertIs(isclose(result[2], -116669/49706, abs_tol=EPSILON), True)


if __name__ == '__main__':
    unittest.main()
