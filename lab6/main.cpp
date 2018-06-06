#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;

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

double abs(double x) {
	return (x > 0) ? x : -x;
}

double norminf(vector <double> a) {
	double mx = 0;
	for (auto &x : a) {
		if (abs(x) > mx) {
			mx = abs(x);
		}
	}
	return mx;
}

int sub(vector <double> &c, vector <double> &a, vector <double> &b) {
	if (a.size() != b.size()) {
		return 1;
	}
	c.resize(a.size());
	for (int i = 0; i < a.size(); i++) {
		c[i] = a[i] - b[i];
	}
	return 0;
}

const double EPS = 1e-4;

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

int main() {
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
	
	vector <double> x0, x1, x2;

	jacobi(x0, a, b);
	gauss(x1, a, b);
	sor(x2, a, b);
	for (auto &x : x0) {
		printf("%.5f ", x);
	}
	printf("\n");
	for (auto &x : x1) {
		printf("%.5f ", x);
	}
	printf("\n");
	for (auto &x : x2) {
		printf("%.5f ", x);
	}
	printf("\n");
}
