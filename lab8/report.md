# 数值分析实验报告

计54 马子轩 2015012283

## 实验要求

用幂法求出下列矩阵按模最大的特征值$\lambda_1$及其对应的特征向量$x_1$，使得$|(\lambda_1)_{k+1} − (\lambda_1)_k|<10^{−5}$
$$
A=\left[
 \begin{matrix}
 5&-4&1\\
 -4&6&-4\\
 1&-4&7
  \end{matrix}
  \right]\tag{1}
$$

$$
B=\left[
 \begin{matrix}
 25&-41&10&-6\\
 -41&68&-17&10\\
 10&-17&5&-3\\
 -6&10&-3&2
  \end{matrix}
  \right]\tag{2}
$$



## 算法描述

幂法的运行过程就是用向量不断迭代乘矩阵的过程，同时为了保持精度，向量每次除以自己的无穷范数，使之归一化，这样不至于几轮迭代之后向量变成0或无穷。

实际实现的时候，先归一化，然后乘矩阵，与之前的特征值进行比较，直到满足条件，输出结果结束。

```cpp
int eigen(vector <double> &vec, double &val, vector <vector <double>> &mat) {
	int n = mat.size();
	for (const auto &row: mat) {
		if (row.size() != n) {
			return 1;
		}
	}
	vector <double> tmp1(n, 1), tmp2(n);
	double delta = 1;
	while (delta > eps) {
		double t = norminf(tmp1);
		for (int i = 0; i < n; i++) {
			tmp2[i] = tmp1[i] / t;
		}
		mul(vec, mat, tmp2);
		delta = abs(norminf(vec) - norminf(tmp1));
		for (int i = 0; i < n; i++) {
			tmp1[i] = vec[i];
		}
	}
	val = norminf(vec);
	for (auto &x : vec) {
		x /= val;
	}
}
```

## 运行结果

```shell
$./main.out < test1.in
0.674020 -1.000000 0.889560
12.254321

$./main.out < test2.in
-0.603972 1.000000 -0.251135 0.148953
98.521698
```

其中test1.in和test2.in是两个输入矩阵.

## 程序清单

main.cpp

```shell
$ g++ main.cpp --std=c++11 -o main.out
$ ./main.out < matrix.in
```



## 体会与展望

本实验代码部分较为简单，而且在之前实验实现的函数基础上，只要对课程内容稍有了解，就能够完成此次实验。
