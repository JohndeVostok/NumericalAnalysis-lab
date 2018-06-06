#include <cstdio>
#include <cstring>
#include <vector>
#include <cmath>

using namespace std;

const double eps = 1e-5;

double tabs(double a) {
	return (a > 0) ? a : -a;
}

double f1(double x) {
	return 2 * x * x * x - x * x + 3 * x - 1;
}

double df1(double x) {
	return 6 * x * x - 2 * x + 3;
}

double f2(double x) {
	return x * x * x - x - 1;
}

double df2(double x) {
	return 3 * x * x - 1;
}

void bs1() {
	double l = -3, r = 3;
	vector <double> res;
	while (r - l > eps) {
		double mid = (l + r) / 2;
		if (f1(mid) <= 0) {
			l = mid;
		} else {
			r = mid;
		}
		res.emplace_back(mid);
	}
	for (const auto &x : res) {
		printf("%.6f ", x);
	}
	printf("%d\n", res.size());
}

void nt1() {
	double x = 0;
	vector <double> res;
	while (tabs(f1(x) / df1(x)) > eps) {
		x -= f1(x) / df1(x);
		res.emplace_back(x);
	}
	for (const auto &x : res) {
		printf("%.6f ", x);
	}
	printf("%d\n", res.size());
}

void test1() {
	bs1();
	nt1();
}

void test21() {
	double x0 = 0, x1 = 0;
	for (int i = 0; i < 10; i++) {
		x0 = pow((x0 + 1) / 2, 1 / 3.0);
		x1 = 2 * x1 * x1 * x1 - 1;
	}
	printf("%.6f, %.6f\n", x0, x1);
}

void nt2(double x0) {
	double x = x0;
	vector <double> res;
	while (tabs(f2(x) / df2(x)) > eps) {
		x -= f2(x) / df2(x);
		res.emplace_back(x);
	}
	for (const auto &x : res) {
		printf("%.6f ", x);
	}
	printf("%d\n", res.size());
}

void test22() {
	nt2(0);
	nt2(1.5);
}

int main() {
	test1();
	test21();
	test22();
}
