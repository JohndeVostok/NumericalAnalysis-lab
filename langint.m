function y = langint(X, Y, x)
    l = size(X, 2);
    y = 0.0;
    for i = 1 : l
        t = 1.0;
        for j = 1 : l
            if i ~= j
                t = t .* (x - X(j)) ./ (X(i) - X(j));
            end
        end
        y = y + t .* Y(i);
    end
end