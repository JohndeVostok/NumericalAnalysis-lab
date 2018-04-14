#include <cstdio>
#include <cstring>
#include <vector>
#include <cmath>

#include "gauss.h"

using namespace std;

int spot(double &c, vector <double> &a, vector <double> &b) {
	if (a.size() != b.size()) {
		return 1;
	}

	c = 0;
	for (int i = 0; i < a.size(); i++) {
		c += a[i] * b[i];
	}
	return 0;
}

double mod(vector <double> &a) {
	double x = 0;
	for (auto i : a) {
		x += i * i;
	}
	return sqrt(x);
}

void init(vector <vector <double>> &x, vector <double> &y) {
	x.clear();
	y.clear();
	x.emplace_back(0);
	for (int i = 0; i < 7; i++) {
		x[0].emplace_back(1);
	}

	x.emplace_back(0);
	x[1].emplace_back(-1.0);
	x[1].emplace_back(-0.5);
	x[1].emplace_back(0.0);
	x[1].emplace_back(0.5);
	x[1].emplace_back(1.0);
	x[1].emplace_back(1.5);
	x[1].emplace_back(2.0);
	
	for (int i = 2; i < 4; i++) {
		x.emplace_back(0);
		for (int j = 0; j < 7; j++) {
			x[i].emplace_back(x[i - 1][j] * x[1][j]);
		}
	}

	y.emplace_back(-4.467);
	y.emplace_back(-0.452);
	y.emplace_back(0.551);
	y.emplace_back(0.048);
	y.emplace_back(-0.447);
	y.emplace_back(0.549);
	y.emplace_back(4.552);
}

void calc(vector <double> &z, vector <vector <double>> &x, vector <double> &y, int dim) {
	vector <vector <double>> ma;
	vector <double> vb;
	double tmp;

	z.clear();
	for (int i = 0; i <= dim; i++) {
		ma.emplace_back(0);
		for (int j = 0; j <= dim; j++) {
			spot(tmp, x[i], x[j]);
			ma[i].emplace_back(tmp);
		}
	}

	for (int i = 0; i <= dim; i++) {
		spot(tmp, x[i], y);
		vb.emplace_back(tmp);
	}

	gauss(z, ma, vb);
}

int main() {
	vector <vector <double>> x;
	vector <double> y;
	vector <double> z;
	init(x, y);
	calc(z, x, y, 1);
	for (auto i : z) {
		printf("%.5f ", i);
	}
	printf("\n");

	calc(z, x, y, 2);
	for (auto i : z) {
		printf("%.5f ", i);
	}
	printf("\n");

	calc(z, x, y, 3);
	for (auto i : z) {
		printf("%.5f ", i);
	}
	printf("\n");
}
