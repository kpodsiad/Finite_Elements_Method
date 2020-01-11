from fem_util import fem

a=0
b=1
f_x = lambda x: 1

n=3
h=(b-a)/n


print(fem(f_x, debug=True))
