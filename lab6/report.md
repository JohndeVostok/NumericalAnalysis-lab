# 数值分析实验报告

计54 马子轩 2015012283

## 实验要求

考虑常微分方程的两点边值问题

$\begin{cases}\epsilon\frac{d^2y}{d^2x}+\frac{dy}{dx}=\alpha\\y(0)=0\\y(1)=1\end{cases},0<a<1$

其精确解为

$y=\frac{1-a}{1-e^{-\frac1\epsilon}}(1-e^{\frac{-x}\epsilon}) + ax$

对$[0,1]$区间n等分

$h=\frac1n,x_i=ih,i=1,2,\dots,n-1$

得到方程

$(\epsilon+h)y_{i+1}-(2\epsilon+h)y_i+\epsilon y_{i-1}=ah^2$

对于$\epsilon=1,a=\frac12,n=100$

用Jacobi,Gauss-Seidel,SOR求解线性方程组,要求四位有效数字,然后和精确解进行比较,并比较上述三种方法.

## 算法描述

### 初始化

迭代法是针对稀疏矩阵比较好的方法,如果按照矩阵进行存储,那么迭代法节省计算量的优点就没有体现,因此我使用图进行矩阵存储.

```cpp
class Node {
public:
	int y, next;
	double s;
	Node (int ty, double ts, int tnext) {
		y = ty;
		s = ts;
		next = tnext;
	}
};

class Matrix {
public:
	vector <Node> nodes;
	vector <int> idx;

	Matrix(int n) {
		idx.resize(n, -1);
	}

	~Matrix() {}

	int addNode(int x, int y, int s) {
		if (x < 0 || x >= idx.size() || y < 0 || y >= idx.size()) {
			return 1;
		}
		nodes.emplace_back(y, s, idx[x]);
		idx[x] = nodes.size() - 1;
		return 0;
	}
};
```

针对题目中矩阵特点进行建图

```cpp
int n = 100;
double eps = 1,	alp = 0.5, h = 1 / double(n);

Matrix a(n - 1);
for (int i = 0; i < n - 1; i++) {
	if (i > 0) {
		a.addNode(i, i - 1, eps);
	}
	a.addNode(i, i, -(2 * eps + h));
	if (i < n) {
		a.addNode(i, i + 1, eps + h);
	}
}
vector <double> b(n - 1, alp * h * h);
b[n - 2] -= eps + h;
```

### Jacobi

Jacobi迭代过程就是使用上一次的结果作为系数直接计算下一轮，直到两次差的无穷范数小于eps.

```cpp
int jacobi(vector <double> &x, Matrix &a, vector <double> &b) {
	if (a.idx.size() != b.size()) {
		return 1;
	}
	vector <double> xp(b.size(), 0), xq(b.size(), 0), dx;
	int gen = 0;
	while (true) {
		gen++;
		for (int i = 0; i < xq.size(); i++) {
			double ci, tmp = b[i];
			for (int j = a.idx[i]; j != -1; j = a.nodes[j].next) {
				if (a.nodes[j].y == i) {
					ci = a.nodes[j].s;
				} else {
					tmp -= xp[a.nodes[j].y] * a.nodes[j].s;
				}
			}
			xq[i] = tmp / ci;
		}
		sub(dx, xp, xq);
		if (norminf(dx) < EPS) {
			break;
		}
		for (int i = 0; i < xp.size(); i++) {
			xp[i] = xq[i];
		}
	}
	printf("%d\n", gen);
	x.clear();
	for (const auto &i : xq) {
		x.emplace_back(i);
	}
	return 0;
}
```

### Gauss-Seidel

Gauss-Seidel方法是针对gauss算法的一个改进，他在计算第k轮的时候使用已经计算出的第k轮结果,和其他的点的k-1轮的结果进行计算.两种方法代码仅有少量不同.

```cpp
int gauss(vector <double> &x, Matrix &a, vector <double> &b) {
	if (a.idx.size() != b.size()) {
		return 1;
	}
	vector <double> xp(b.size(), 0), xq(b.size(), 0), dx;
	int gen = 0;
	while (true) {
		gen++;
		for (int i = 0; i < xq.size(); i++) {
			double ci, tmp = b[i];
			for (int j = a.idx[i]; j != -1; j = a.nodes[j].next) {
				if (a.nodes[j].y == i) {
					ci = a.nodes[j].s;
				} else {
					if (a.nodes[j].y < i) {
						tmp -= xq[a.nodes[j].y] * a.nodes[j].s;
					} else {
						tmp -= xp[a.nodes[j].y] * a.nodes[j].s;
					}
				}
			}
			xq[i] = tmp / ci;
		}
		sub(dx, xp, xq);
		if (norminf(dx) < EPS) {
			break;
		}
		for (int i = 0; i < xp.size(); i++) {
			xp[i] = xq[i];
		}
	}
	printf("%d\n", gen);
	x.clear();
	for (const auto &i : xq) {
		x.emplace_back(i);
	}
	return 0;
}
```

### SOR

SOR方法是针对Gauss-Seidel方法的一种改进,引进了松弛系数$\omega$,松弛系数可以理解为变化的权重.在Gauss-Seidel算法中,两次x的差乘松弛系数,加到上一轮的结果上,得到下一轮.这样,通过调整松弛系数,即可调整Gauss-Seidel的权重,显然松弛系数为1的时候,SOR和Gauss-Seidel算法一致,松弛系数小于1的时候.Gauss-Seidel的调整能力降低,而松弛系数大于1的时候,Gauss-Seidel的调整能力高.这种情况下,通过选取合适的$\omega$可以在更少的轮数中迭代出结果.

在实验过程中,经过测试在w取1.933的时候,迭代186轮即可得到结果.优化效果显著.

```cpp
int sor(vector <double> &x, Matrix &a, vector <double> &b) {
	if (a.idx.size() != b.size()) {
		return 1;
	}
	double w = 1.933;
	vector <double> xp(b.size(), 0), xq(b.size(), 0), dx;
	int gen = 0;
	while (true) {
		gen++;
		for (int i = 0; i < xq.size(); i++) {
			double ci, tmp = b[i];
			for (int j = a.idx[i]; j != -1; j = a.nodes[j].next) {
				if (a.nodes[j].y == i) {
					ci = a.nodes[j].s;
				} else {
					if (a.nodes[j].y < i) {
						tmp -= xq[a.nodes[j].y] * a.nodes[j].s;
					} else {
						tmp -= xp[a.nodes[j].y] * a.nodes[j].s;
					}
				}
			}
			xq[i] = (1 - w) * xp[i] + w * tmp / ci;
		}
		sub(dx, xp, xq);
		if (norminf(dx) < EPS) {
			break;
		}
		for (int i = 0; i < xp.size(); i++) {
			xp[i] = xq[i];
		}
	}
	printf("%d\n", gen);
	x.clear();
	for (const auto &i : xq) {
		x.emplace_back(i);
	}
	return 0;
}
```

## 运行结果

迭代算法的运行结果如下,从上到下依次是Jacobi的迭代次数,Gauss-Seidel的迭代次数,SOR的迭代次数.Jacobi的结果,Gauss-Seidel的结果,SOR的结果.

```shell
3641
1793
186
0.00463 0.00931 0.01405 0.01884 0.02370 0.02862 0.03361 0.03866 0.04381 0.04900 0.05430 0.05965 0.06513 0.07064 0.07630 0.08199 0.08784 0.09372 0.09976 0.10585 0.11210 0.11839 0.12487 0.13138 0.13808 0.14482 0.15175 0.15873 0.16591 0.17313 0.18056 0.18803 0.19572 0.20345 0.21140 0.21939 0.22762 0.23588 0.24438 0.25292 0.26169 0.27051 0.27957 0.28867 0.29801 0.30739 0.31703 0.32670 0.33662 0.34658 0.35679 0.36704 0.37754 0.38807 0.39886 0.40969 0.42076 0.43187 0.44323 0.45463 0.46627 0.47795 0.48986 0.50182 0.51401 0.52624 0.53869 0.55119 0.56391 0.57667 0.58964 0.60265 0.61587 0.62913 0.64259 0.65610 0.66979 0.68352 0.69743 0.71139 0.72551 0.73968 0.75401 0.76838 0.78291 0.79747 0.81217 0.82692 0.84180 0.85672 0.87175 0.88683 0.90201 0.91724 0.93255 0.94792 0.96336 0.97885 0.99440 
0.00439 0.00883 0.01333 0.01790 0.02253 0.02723 0.03200 0.03685 0.04177 0.04678 0.05187 0.05704 0.06231 0.06767 0.07312 0.07867 0.08432 0.09008 0.09594 0.10190 0.10798 0.11418 0.12049 0.12691 0.13346 0.14013 0.14692 0.15385 0.16090 0.16808 0.17539 0.18284 0.19042 0.19814 0.20600 0.21400 0.22214 0.23042 0.23885 0.24742 0.25613 0.26500 0.27401 0.28317 0.29248 0.30193 0.31154 0.32130 0.33120 0.34126 0.35147 0.36182 0.37233 0.38298 0.39379 0.40474 0.41584 0.42709 0.43848 0.45002 0.46170 0.47353 0.48550 0.49760 0.50985 0.52224 0.53476 0.54742 0.56021 0.57313 0.58618 0.59936 0.61266 0.62609 0.63963 0.65330 0.66708 0.68098 0.69498 0.70910 0.72332 0.73765 0.75208 0.76660 0.78123 0.79594 0.81074 0.82564 0.84061 0.85567 0.87080 0.88601 0.90129 0.91664 0.93205 0.94753 0.96307 0.97866 0.99430 
0.00749 0.01504 0.02264 0.03030 0.03801 0.04578 0.05359 0.06147 0.06939 0.07737 0.08540 0.09349 0.10162 0.10987 0.11813 0.12645 0.13482 0.14324 0.15172 0.16025 0.16884 0.17748 0.18618 0.19492 0.20373 0.21258 0.22149 0.23045 0.23947 0.24854 0.25766 0.26683 0.27606 0.28534 0.29467 0.30406 0.31349 0.32298 0.33252 0.34212 0.35176 0.36146 0.37120 0.38100 0.39086 0.40076 0.41071 0.42072 0.43077 0.44088 0.45104 0.46125 0.47151 0.48182 0.49218 0.50259 0.51305 0.52357 0.53413 0.54474 0.55540 0.56612 0.57688 0.58769 0.59856 0.60947 0.62043 0.63145 0.64251 0.65362 0.66479 0.67600 0.68726 0.69857 0.70993 0.72134 0.73280 0.74431 0.75587 0.76748 0.77913 0.79084 0.80260 0.81440 0.82626 0.83816 0.85012 0.86212 0.87417 0.88627 0.89842 0.91063 0.92287 0.93517 0.94752 0.95992 0.97237 0.98486 0.99741 
```

我们可以看出,SOR方法在松弛系数选择得当的情况下有明显的性能优势.

我们针对结果进行分析

```python
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
```

运行结果

```shell
0.2148378972172078
0.22023789721720782
0.12090883835280009
```

可以看出,前两种方法误差相近,SOR方法的误差更小.

考虑到误差来源,一方面来源于迭代法的误差,另一方面来源于离散化的误差.而离散化部分的相对误差达到了$10^{-2}$数量级,因此最终误差的无穷范数在这个数量是可以理解的.

综合考虑,在这个问题上Gauss-Seidel方法优于Jacobi方法.而SOR方法在松弛系数选取得当的情况下明显优于以上两种方法.但是$\omega$的选取需要额外进行较多的实验,实验中取小于1的数值时候,要劣于Gauss-Seidel方法.而取到1.5以上之前性能优化都没有这么明显.因此,在真实使用的情况下,SOR方法不一定能够节省时间,提高效率.

## 程序清单

main.cpp

res.py

```shell
g++ main.cpp --std=c++11 -o main
./main
./main > res.out
python res.py
```

其中main.cpp为计算部分

res.py为误差分析

## 体会与展望

本次实验中结果与我预期相差较大,其中在离散化误差上花费了比较多的时间进行思考求证.才得出这一部分误差影响非常大的结论.但是针对计算部分,我认为这种处理稀疏矩阵的方式,是非常合理的.能够作为后续实验的理论指导.