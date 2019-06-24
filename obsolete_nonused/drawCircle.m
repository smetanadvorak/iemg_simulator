function [x_circle, y_circle] = drawCircle(center, rad, ax)

x_circle = linspace(-rad, rad, 100);
y_circle = sqrt(rad^2 - x_circle.^2);
x_circle = [x_circle, x_circle(end:-1:1)];
y_circle = [y_circle, -y_circle(end:-1:1)];

x_circle = x_circle + center(1);
y_circle = y_circle + center(2);

if nargin >= 3
    axes(ax);
    plot(x_circle, y_circle, 'b');
end