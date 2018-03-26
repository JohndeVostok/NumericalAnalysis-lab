#include <cstdio>

const int n0 = 10, n1 = 20;
const float min = -5, max = 5;

float x0[11], x1[21], y0[11], y1[21];

float f(float x) {
	return 1 / (1 + 16 * x * x);
}

float l0(float x) {
	float y = 0, t = 1;
	for (int i = 0; i <= n0; i++) {
		t = 1;
		for (int j = 0; j <= n0; j++) {
			if (i != j) {
				t = t * (x - x0[j])/(x0[i] - x0[j]);
			}
		}
		y = y + t * y[i];
	}
	return y;
}

float l1(float x) {
	float y = 0, t = 1;
	for (int i = 0; i <= n1; i++) {
		t = 1;
		for (int j = 0; j <= n1; j++) {
			if (i != j) {
				t = t * (x - x1[j])/(x1[i] - x1[j]);
			}
		}
		y = y + t * y[i];
	}
	return y;
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
}
