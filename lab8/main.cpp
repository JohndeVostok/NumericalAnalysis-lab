#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;

const double eps = 1e-5;
const double inf = 1e8;

double abs(double a) {
	return (a > 0) ? a : -a;
}

int mul(vector <double> &res, vector <vector <double>> &mat, vector <double> &vec) {
	int n = vec.size();
	//validate matrix and vector
	if (mat.size() != n) {
		return 1;
	}
	for (const auto &row : mat) {
		if (row.size() != n) {
			return 1;
		}
	}
	res.resize(n);
	for (int i = 0; i < n; i++) {
		res[i] = 0;
		for (int j = 0; j < n; j++) {
			res[i] += mat[i][j] * vec[j];
		}
	}
}

double norminf(vector <double> &vec) {
	double res = -inf;
	for (const auto &x : vec) {
		if (abs(x) > res) {
			res = abs(x);
		}
	}
	return res;
}

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

int main() {
	int n;
	vector <vector <double>> mat;
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		mat.emplace_back(0);
		for (int j = 0; j < n; j++) {
			double x;
			scanf("%lf", &x);
			mat[i].emplace_back(x);
		}
	}
	vector <double> vec;
	double val;
	eigen(vec, val, mat);
	for (const auto &x : vec) {
		printf("%f ", x);
	}
	printf("\n%f\n", val);
}
