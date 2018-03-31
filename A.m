X = [-5: 0.001: 5];
Y = 1 ./ (1 + 16 * X .^ 2);

x1 = [-5: 0.1: 5];
x2 = [-5: 0.05: 5];

y1 = 1 ./ (1 + 16 * x1 .^ 2);
y2 = 1 ./ (1 + 16 * x2 .^ 2);

LY1 = langint(x1, y1, X);
LY2 = langint(x2, y2, X);

SY1 = spline(x1, y1, X);
SY2 = spline(x2, y2, X);

plot(X, Y);
hold on;
plot(X, SY1);
hold on;
plot(X, SY2);
