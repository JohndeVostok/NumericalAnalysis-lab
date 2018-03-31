X = [-5: 0.01: 5];
Y = 1 ./ (1 + 16 * X .^ 2);

x1 = [-5: 1: 5];
x2 = [-5: 0.5: 5];

y1 = 1 ./ (1 + 16 * x1 .^ 2);
y2 = 1 ./ (1 + 16 * x2 .^ 2);

LY1 = langint(x1, y1, X);
LY2 = langint(x2, y2, X);

plot(X, Y);
hold on;
plot(X, LY1);
plot(X, LY2);
