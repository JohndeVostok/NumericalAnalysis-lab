#include <cstdio>
#include <cstring>

const int n0 = 10, n1 = 20;
const double min = -5, max = 5, eps = 1e-5;

int  c[80], d[80];
double x0[11], x1[21], y0[11], y1[21];
double a[80][80], b[80], e[80];
double s0[10][4], s1[20][4];

double abs(double x) {
	return (x > 0) ? x : -x; 
}

double f(double x) {
	return 1 / (1 + 16 * x * x);
}

double l0(double x) {
	double ans = 0, t = 1;
	for (int i = 0; i <= n0; i++) {
		t = 1;
		for (int j = 0; j <= n0; j++) {
			if (i != j) {
				t = t * (x - x0[j])/(x0[i] - x0[j]);
			}
		}
		ans = ans + t * y0[i];
	}
	return ans;
}

double l1(double x) {
	double ans = 0, t = 1;
	for (int i = 0; i <= n1; i++) {
		t = 1;
		for (int j = 0; j <= n1; j++) {
			if (i != j) {
				t = t * (x - x1[j])/(x1[i] - x1[j]);
			}
		}
		ans = ans + t * y1[i];
	}
	return ans;
}

void print(int n) {
	printf("erased:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%.2f ", a[i][j]);
		}
		printf("%f\n", b[i]);
	}
	printf("\n");
}

bool gauss(int n) {
	int p;
	double t, max;
	for (int i = 0; i < n; i++) {
		p = -1;
		max = 0;
		for (int j = 0; j < n; j++) {
			if (!c[j] && abs(a[j][i]) > max) {
				p = j;
				max = abs(a[j][i]);
			}
		}
		if (p == -1) return 1;
		t = a[p][i];
		for (int j = 0; j < n; j++) {
			a[p][j] /= t;
		}
		b[p] /= t;
		c[p] = 1;
		d[i] = p;
		for (int j = 0; j < n; j++) {
			if (p != j) {
				t = a[j][i];
				for (int k = 0; k < n; k++) {
					a[j][k] -= t * a[p][k];
				}
				b[j] -= t * b[p];
			}
		}
	}
	for (int i = 0; i < n; i++) {
		e[i] = b[d[i]] / a[d[i]][i];
	}
	return 0;
}

int main() {
	for (int i = 0; i <= n0; i++) {
		x0[i] = min + (max - min) * i / n0;
		y0[i] = f(x0[i]);
	}

	for (int i = 0; i < n1; i++) {
		x1[i] = min + (max - min) * i / n1;
		y1[i] = f(x1[i]);
	}
	
	//S0
	memset(a, 0, sizeof(a));
	memset(b, 0, sizeof(b));
	memset(c, 0, sizeof(c));
	memset(d, 0, sizeof(d));
	memset(e, 0, sizeof(e));
	int cnt = 0;
	for (int i = 0; i < n0; i++) {
		a[cnt][i * 4] = 1;
		a[cnt][i * 4 + 1] = x0[i];
		a[cnt][i * 4 + 2] = x0[i] * x0[i];
		a[cnt][i * 4 + 3] = x0[i] * x0[i] * x0[i];
		b[cnt] = y0[i];
		cnt++;
		a[cnt][i * 4] = 1;
		a[cnt][i * 4 + 1] = x0[i + 1];
		a[cnt][i * 4 + 2] = x0[i + 1] * x0[i + 1];
		a[cnt][i * 4 + 3] = x0[i + 1] * x0[i + 1] * x0[i + 1];
		b[cnt] = y0[i + 1];
		cnt++;
	}
	for (int i = 1; i < n0; i++) {
		a[cnt][i * 4 + 1] = 1;
		a[cnt][i * 4 + 2] = 2 * x0[i];
		a[cnt][i * 4 + 3] = 3 * x0[i] * x0[i];
		a[cnt][i * 4 - 3] = -1;
		a[cnt][i * 4 - 2] = -2 * x0[i];
		a[cnt][i * 4 - 1] = -3 * x0[i] * x0[i];
		cnt++;
	}
	for (int i = 1; i < n0; i++) {
		a[cnt][i * 4 + 2] = 2;
		a[cnt][i * 4 + 3] = 6 * x0[i];
		a[cnt][i * 4 - 2] = -2;
		a[cnt][i * 4 - 1] = -6 * x0[i];
		cnt++;
	}
	a[cnt][1] = 1;
	a[cnt][2] = 2 * x0[0];
	a[cnt][3] = 3 * x0[0] * x0[0];
	b[cnt] = 0.000995;
	cnt++;
	a[cnt][4 * n0 - 3] = 1;
	a[cnt][4 * n0 - 2] = 2 * x0[n0];
	a[cnt][4 * n0 - 1] = 3 * x0[n0] * x0[n0];
	b[cnt] = -0.000995;
	cnt++;

	FILE *f = fopen("s0.txt", "w");
	bool err = gauss(40);
	for (int i = 0; i < n0; i++) {
		for (int j = 0; j < 4; j++) {
			s0[i][j] = e[i * 4 + j];
			fprintf(f, "%f ", s0[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	
	//S1
	memset(a, 0, sizeof(a));
	memset(b, 0, sizeof(b));
	memset(c, 0, sizeof(c));
	memset(d, 0, sizeof(d));
	memset(e, 0, sizeof(e));
	cnt = 0;
	for (int i = 0; i < n1; i++) {
		a[cnt][i * 4] = 1;
		a[cnt][i * 4 + 1] = x1[i];
		a[cnt][i * 4 + 2] = x1[i] * x1[i];
		a[cnt][i * 4 + 3] = x1[i] * x1[i] * x1[i];
		b[cnt] = y1[i];
		cnt++;
		a[cnt][i * 4] = 1;
		a[cnt][i * 4 + 1] = x1[i + 1];
		a[cnt][i * 4 + 2] = x1[i + 1] * x1[i + 1];
		a[cnt][i * 4 + 3] = x1[i + 1] * x1[i + 1] * x1[i + 1];
		b[cnt] = y1[i + 1];
		cnt++;
	}
	for (int i = 1; i < n1; i++) {
		a[cnt][i * 4 + 1] = 1;
		a[cnt][i * 4 + 2] = 2 * x1[i];
		a[cnt][i * 4 + 3] = 3 * x1[i] * x1[i];
		a[cnt][i * 4 - 3] = -1;
		a[cnt][i * 4 - 2] = -2 * x1[i];
		a[cnt][i * 4 - 1] = -3 * x1[i] * x1[i];
		cnt++;
	}
	for (int i = 1; i < n1; i++) {
		a[cnt][i * 4 + 2] = 2;
		a[cnt][i * 4 + 3] = 6 * x1[i];
		a[cnt][i * 4 - 2] = -2;
		a[cnt][i * 4 - 1] = -6 * x1[i];
		cnt++;
	}
	a[cnt][1] = 1;
	a[cnt][2] = 2 * x1[0];
	a[cnt][3] = 3 * x1[0] * x1[0];
	b[cnt] = 0.000995;
	cnt++;
	a[cnt][4 * n1 - 3] = 1;
	a[cnt][4 * n1 - 2] = 2 * x1[n0];
	a[cnt][4 * n1 - 1] = 3 * x1[n0] * x1[n1];
	b[cnt] = -0.000995;
	cnt++;

	f = fopen("s1.txt", "w");
	err = gauss(4 * n1);
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < 4; j++) {
			s1[i][j] = e[i * 4 + j];
			fprintf(f, "%f ", s1[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

