#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

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

int main() {
	double ans1 = calc1();
	printf("trapezoid: %.8f\n", ans1);
	double ans2 = calc2();
	printf("Simpson: %.8f\n", ans2);
}
