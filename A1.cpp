#include <cstdio>

const float ln2 = 0.693147190546;

float abs(float x) {
	return (x > 0) ? x : -x;
}

int main() {
	int n;
	float eps = 0.00005, x = 0.0;
	for (n = 1; abs(x - ln2) > eps; x += (n & 1) ? 1.0 / n : -1.0 / n, n++);
	printf("%d %.8f\n", n, x);
}
