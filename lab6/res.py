import math

eps = 1
a = 0.5
n = 100
h = 1 / n

def f(x):
	y = (1 - a) / (1 - math.exp(-1 / eps)) * (1 - math.exp(-x / eps)) + a * x;
	return y

if __name__ == "__main__":
	y = [f((i + 1) / n) for i in range(99)]
	with open("res.out", "r") as f:
		lines = f.readlines()
	yj = [float(x) for x in lines[3].split(" ")[:-1]]
	yg = [float(x) for x in lines[4].split(" ")[:-1]]
	ys = [float(x) for x in lines[5].split(" ")[:-1]]

	dj = max([abs(y[i] - yj[i]) for i in range(99)])
	print(dj)
	dg = max([abs(y[i] - yg[i]) for i in range(99)])
	print(dg)
	ds = max([abs(y[i] - ys[i]) for i in range(99)])
	print(ds)
