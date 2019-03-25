import random

def witness(a, n):
    t, u = 0, n - 1
    while u % 2 == 0:
        t += 1
        u >>= 1
    x, _x = pow(a, u, n), 0
    while t:
        _x = (x * x) % n
        if _x == 1 and x != 1 and x != n - 1:
            return True
        x = _x
        t -= 1
    return _x != 1

def isPrime(x):
    if x < 2:
        return False
    for i in {2, 3, 5, 7, 11, 13, 17, 19}:
        if x % i == 0:
            return x == i
    for i in range(10):
        if witness(random.randint(1, x), x):
            return False
    return True


for i in range(1, 1000):
    if isPrime((2 ** i) + 1):
        print(i, (2 ** i) + 1)
