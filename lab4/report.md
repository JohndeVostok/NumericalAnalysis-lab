# 数值分析实验报告

计54 马子轩 2015012283

## 实验内容

1. 用不同数值方法计算积分$I=\int_0^1e^xdx$

(1) 使用复合梯形公式，并回答分为多少份可以使误差不超过$1\times10^{-6}$.

(2) 使用复合辛普森公式，并回答分为多少份可以达到相同精度.

(3) 使用龙贝格积分，并回答分为多少份可以达到相同精度.

2. 用下面的复合高斯公式做近似积分$\pi=4\int_0^1\frac{1}{1+x^2}dx$.

将$[a,b]$做等距划分,$x_i=a+ih(i=0,...,n), h=(b-a)/n$.

在每个区间内用两点高斯公式，有以下结果

$\int_a^bf(x)dx=\frac{h}{2}\sum_{i=0}^{n-1}[f(x_{i+\frac{1}{2}}-\frac{h}{2\sqrt{3}})+f(x_{i+\frac{1}{2}}+\frac{h}{2\sqrt{3}})]+\frac{(b-a)h^4}{4320}f^{(4)}(\xi),(\xi\in(a,b))$

其中$x_{i+\frac{1}{2}}=x_i+\frac{h}{2}​$,对h做先验估计，再用上式近似积分,将理论值与实验结果比较.

## 实验过程

### 复合梯形公式

$R_n(f)=-\frac{b-a}{12}h^2f''(\eta),\eta\in(a,b)$

$f''(\eta)=(e^\eta)''=e^\eta<e^1=e

$1\times10^{-6}>|R_n(f)|=\frac{1}{12}h^2f''(\eta)$

$1\times10^{-6} >\frac{1}{12}h^2e$

$h<\sqrt{12\times1\times10^{-6}\div e}=0.0021011$

$n\ge\frac{1}{h}=475.94,n=476$

实验代码

```cpp
double calc1() {
	int n = 476;
	double h = 1 / double(n);
	vector <double> x, y;
	for (int i = 0; i <= n; i++) {
		x.emplace_back(i / double(n));
		y.emplace_back(exp(x.back()));
	}
	double res = 0;
	for (int i = 0; i < n; i++) {
		res += h / 2 * (y[i + 1] + y[i]);
	}
	return res;
}
```

运行结果

```shell
trapezoid: 1.71828246
```

误差小于$1\times10^{-6}$

### 复合辛普森公式

$R_n(f)=-\frac{b-a}{180}(\frac{h}{2})^4f^{(4)}(\eta),\eta\in(a,b)$

$f^{(4)}(\eta)=e^\eta<e$

$1\times10^{-6}>|R_n(f)|=\frac{1}{2880}h^4f^{(4)}(\eta)$

$1\times10^{-6}>\frac{1}{2880}h^4e$

$h<(2880 \times1\times10^{-6}\div e)^{\frac{1}{4}}=0.2316584$

$n\ge\frac{1}{h},n=5$

实验代码

```cpp
double calc2() {
	int n = 5;
	double h = 1 / double(n);
	vector <double> x, y;
	for (int i = 0; i <= 2 * n; i++) {
		x.emplace_back(i / double(2 * n));
		y.emplace_back(exp(x.back()));
	}
	double res = 0;
	for (int i = 0; i < n; i++) {
		res += h / 6 * (y[2 * i] + 4 * y[2 * i + 1] + y[2 * i + 2]);
	}
	return res;
}
```

```shell
Simpson: 1.71828278
```

误差小于$1\times10^{-6}$

### 龙贝格求积公式

实验代码

```cpp
pair <int, double> calc3() {
	int k = 0, n = 1;
	double h = 1, tmp = 1;
	vector <vector <double>> t;
	t.emplace_back(0);
	t[0].emplace_back(0.5 * (exp(1) - 1));
	while (true) {
		k++;
		n = (1 << k);
		h = 1 / double(n);
		t.emplace_back(0);
		t[k].emplace_back(t[k - 1][0] / 2);
		for (int i = 0; i < n; i++) {
			t[k][0] += h / 2 * exp((2 * i + 1) / double(2 * n));
		}
		tmp = 1;
		for (int i = 1; i <= k; i++) {
			tmp *= 4;
			t[k].emplace_back((tmp * t[k][i - 1] - t[k - 1][i - 1]) / (tmp - 1));
		}
		if (abs(t[k][k] - t[k - 1][k - 1]) < eps) {
			break;
		}
	}
	return make_pair(k, t[k][k]);
}
```

```shell
Romberg: 20 1.71828131
```

误差小于$1\times10^{-6}$

迭代次数20，分段数1048576

### 复合高斯公式

$(\frac{1}{1 + x^2})^{(4)}\le4!=24$

$1\times10^{-6}\ge|R_n(f)|=\frac{h^4}{4320}f^{(4)}(\xi)$

$h^4<180\times10^{-6}$

$h<0.1158292$

$n\ge\frac{1}{h}=8.6334002,n=9$

$\int_a^bf(x)dx=0.78539816$

```cpp
double calc4() {
	int n = 9;
	double h = 1 / double(n), res = 0, tmp = 2 * sqrt(3);
	for (int i = 0; i < n; i++) {
		res += h / 2 * f4((2 * i + 1) / double(2 * n) - h / tmp);
		res += h / 2 * f4((2 * i + 1) / double(2 * n) + h / tmp);
	}
	return res;
}
```

```shell
0.78539816
```

误差小于$1\times10^{-6}$

##实验组织

代码清单

main.cpp

运行方法

```shell
g++ main.cpp --std=c++11 -o main
./main
```


## 实验总结

这个实验做起来较为轻松，主要原因是这次不需要画图了。而实际写的时候，也体会到了不同的积分方式效率区别，我比较不理解的是龙贝格积分为什么会比直接梯形还慢，这一点还需要进行考虑。
