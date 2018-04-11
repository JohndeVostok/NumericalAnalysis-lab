#include <cstdio>
#include <cstring>
#include <vector>
#include <cmath>

using namespace std;

vector <double> x, y;

void init() {
	x.emplace_back(-1.0);
	y.emplace_back(-4.467);
	x.emplace_back(-0.5);
	y.emplace_back(-0.452);
	x.emplace_back(0.0);
	y.emplace_back(0.551);
	x.emplace_back(0.5);
	y.emplace_back(0.048);
	x.emplace_back(1.0);
	y.emplace_back(-0.477);
	x.emplace_back(1.5);
	y.emplace_back(0.549);
	x.emplace_back(2.0);
	y.emplace_back(4.552);
}

int spot(vector <double> &c, vector <double> &a, vector <double> &b) {
	if (a.size() != b.size()) {
		return 1;
	}

	c.clear();
	for (int i = 0; i < a.size(); i++) {
		c.emplace_back(a[i] * b[i]);
	}
}

double mod(vector <double> &a) {
	double x = 0;
	for (auto i : a) {
		x += i * i;
	}
	return sqrt(x);
}

void calc1(vector <double> &x, vector <double> &y, double a, double b) {
	
}

int main() {
	double a, b, c, d;
//f(x) = a + bx
	calc1(x, y, a, b);
	printf("f(x) = a + bx :\n");
	printf("%.5fx + %.5f\n", a, b);
}
