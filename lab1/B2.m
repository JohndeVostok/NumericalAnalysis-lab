X = [-5: 0.01: 5];
Y = 1 ./ (1 + 16 * X .^ 2);

x1 = [-5: 1: 5];
x2 = [-5: 0.5: 5];

y1 = 1 ./ (1 + 16 * x1 .^ 2);
y2 = 1 ./ (1 + 16 * x2 .^ 2);

s1 = load("s0.txt");
s2 = load("s1.txt");

for i = 1 : 10
    x = X((i - 1) * 100 + 1 : i * 100);
    sy1((i - 1) * 100 + 1 : i * 100) = s1(i, 4) * x .^ 3 + s1(i, 3) * x .^ 2 + s1(i, 2) * x + s1(i, 1);
end

for i = 1 : 20
    x = X((i - 1) * 50 + 1 : i * 50);
    sy2((i - 1) * 50 + 1 : i * 50) = s2(i, 4) * x .^ 3 + s2(i, 3) * x .^ 2 + s2(i, 2) * x + s2(i, 1);
end

sy1(1001) = Y(1001);
sy2(1001) = Y(1001);

plot(X, Y);
hold on;
plot(X, sy1);
plot(X, sy2);