x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
y = [-4.467, -0.452, 0.551, 0.048, -0.447, 0.549, 4.552];

p1 = [2.00436, -0.95446]
p2 = [-0.00376, 2.00812, -0.95164]
p3 = [2.00356, -3.00910, 0.00456, 0.55102]

plot(x, y, 'o');

hold on;

tx = [-1.0 : 0.01 : 2.0];

plot(tx, polyval(p1, tx));
hold on;
plot(tx, polyval(p2, tx));
hold on;
plot(tx, polyval(p3, tx));

saveas(gcf, 'plot.jpg', 'jpg');