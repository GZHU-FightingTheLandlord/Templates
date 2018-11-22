# 非比赛可以使用numpy时

from scipy.interpolate import lagrange
x = [1, 2, 4, 5]
y = [2, 3, 5, 6]
a = lagrange(x, y)
print(a)

"""
这指数表达方式有点秀
output:
          3             2
1.11e-16 x - 8.882e-16 x + 1 x + 1
"""
