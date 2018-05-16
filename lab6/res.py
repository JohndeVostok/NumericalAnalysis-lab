import math

eps = 1
a = 0.5
n = 100
h = 1 / n

def f(x):
	y = (1 - a) / (1 - math.exp(-1 / eps)) * (1 - math.exp(-x / eps)) + a * x;
	return y

if __name__ == "__main__":
	y = [f((i) / n) for i in range(101)]
	print(y)

	p = eps
	q = -(2 * eps + h)
	r = eps + h

	print(p * y[98] + q * y[99], a * h * h - r)
