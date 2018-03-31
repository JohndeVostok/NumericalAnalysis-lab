#include <cstdio>
#include <cstring>

const int n0 = 10, n1 = 20;
const float min = -5, max = 5, eps = 1e-7;

int d[80];
float x0[11], x1[21], y0[11], y1[21];
float a[80][80], b[80], c[80], e[80];
float s0[10][4], s1[20][4];

float abs(float x) {
	return (x > 0) ? x : -x; 
}

float f(float x) {
	return 1 / (1 + 16 * x * x);
}

float l0(float x) {
	float ans = 0, t = 1;
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

float l1(float x) {
	float ans = 0, t = 1;
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
			printf("%f ", a[i][j]);
		}
		printf("%f\n", b[i]);
	}
	printf("\n");
}

void gauss(int n) {
	int p;
	float t;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (!c[j] && abs(a[j][i]) > eps) {
				p = j;
				break;
			}
		}
		c[p] = 1;
		for (int j = 0; j < n; j++) {
			if (p != j && abs(a[j][i]) > eps) {
				t = a[j][i] / a[p][i];
				for (int k = 0; k < n; k++) {
					a[j][k] -= t * a[p][k];
				}
				b[j] -= t * b[p];
			}
		}
		d[i] = p;
	}
	for (int i = 0; i < n; i++) {
		e[i] = b[d[i]] / a[d[i]][i];
	}
}

int main() {
	for (int i = 0; i <= n0; i++) {
		x0[i] = min + (max - min) * i;
		y0[i] = f(x0[i]);
	}

	for (int i = 0; i < n1; i++) {
		x1[i] = min + (max - min) * i;
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
	gauss(10);
	for (int i = 0; i < n0; i++) {
		for (int j = 0; j < 4; j++) {
			s0[i][j] = e[i * n0 + j];
			printf("%f ", s0[i][j]);
		}
		printf("\n");
	}
}

