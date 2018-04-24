#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

const double eps = 1e-6;

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

double f4(double x) {
	return 1 / (1 + x * x);
}

double calc4() {
	int n = 9;
	double h = 1 / double(n), res = 0, tmp = 2 * sqrt(3);
	for (int i = 0; i < n; i++) {
		res += h / 2 * f4((2 * i + 1) / double(2 * n) - h / tmp);
		res += h / 2 * f4((2 * i + 1) / double(2 * n) + h / tmp);
	}
	return res;
}

int main() {
	double ans1 = calc1();
	printf("trapezoid: %.8f\n", ans1);
	double ans2 = calc2();
	printf("Simpson: %.8f\n", ans2);
	pair <int, double> ans3 = calc3();
	printf("Romberg: %d %.8f\n", ans3.first, ans3.second);
	double ans4 = calc4();
	printf("Gauss: %.8f\n", ans4);
}
