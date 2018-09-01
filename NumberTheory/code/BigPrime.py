import random


def fpow(a, x, n):
	ans = 1
	while x > 0:
		if x & 1:
			ans = ans * a % n
		a = a * a % n
		x >>= 1
	return ans


def check(a, n, x, t):
	ret = fpow(a, x, n)
	last = ret
	for i in range(0, t):
		ret = ret * ret % n
		if ret == 1 and last != 1 and last != n - 1:
			return True
		last = ret
	if ret != 1:
		return True
	return False


# there change the times of Rabin-Miller
TIMES = 50


def is_prime(n):
	if n == 0 or n == 1:
		return False
	if n <= (1 << 20):
		i = 2
		while i * i <= n:
			if n % i == 0:
				return False
			i += 1
		return True
	else:
		if ~n & 1:
			return False
		x = n - 1
		t = 0
		while not x & 1:
			x >>= 1
			t += 1
		for i in range(0, TIMES):
			a = random.randint(1, n - 2)
			if check(a, n, x, t):
				return False
		return True


"""
from BigPrime import is_prime
is_prime(n) return whether n is a prime
"""
if __name__ == '__main__':
	print('Yes' if is_prime(int(input())) else 'No')
