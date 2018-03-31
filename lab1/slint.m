function s = slint(X, Y)
    l = size(X, 2) - 1;
    cnt = 1;
    a = zeros(4 * l, 4 * l);
    b = zeros(4 * l);
	for i = 1 : l
        a(cnt, i * 4 - 3) = 1;
        a(cnt, i * 4 - 2) = X(i);
        a(cnt, i * 4 - 1) = X(i) * X(i);
        a(cnt, i * 4) = X(i) * X(i) * X(i);
        b(cnt) = Y(i);
        cnt = cnt + 1;
        a(cnt, i * 4 - 3) = 1;
        a(cnt, i * 4 - 2) = X(i + 1);
        a(cnt, i * 4 - 1) = X(i + 1) * X(i + 1);
        a(cnt, i * 4) = X(i + 1) * X(i + 1) * X(i + 1);
        b(cnt) = Y(i + 1);
        cnt = cnt + 1;
    end
	for i = 2 : l
		a(cnt, i * 4 - 2) = 1;
		a(cnt, i * 4 - 1) = 2 * X(i);
		a(cnt, i * 4) = 3 * X(i) * X(i);
		a(cnt, i * 4 - 6) = -1;
		a(cnt, i * 4 - 5) = -2 * X(i);
		a(cnt, i * 4 - 4) = -3 * X(i) * X(i);
		cnt = cnt + 1;
    end
    for i = 2 : l
		a(cnt, i * 4 - 1) = 2;
		a(cnt, i * 4) = 6 * X(i);
		a(cnt, i * 4 - 5) = -2;
		a(cnt, i * 4 - 4) = -6 * X(i);
		cnt = cnt + 1;
    end
	a(cnt, 2) = 1;
	a(cnt, 3) = 2 * X(1);
	a(cnt, 4) = 3 * X(1) * X(1);
	b(cnt) = 0.000995;
	cnt = cnt + 1;
	a(cnt, 4 * l - 2) = 1;
	a(cnt, 4 * l - 1) = 2 * X(l + 1);
	a(cnt, 4 * l) = 3 * X(l + 1) * X(l + 1);
	b(cnt) = -0.000995;
    c = inv(a) * b;
    s = zeros(l, 4);
    for i = 1 : l
        for j = 1 : 4
            s(i, j) = c((i - 1) * 4 + j);
        end
    end
end