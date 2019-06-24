% PDF
mu = [-0.5,0];
sigma = 0.1 * eye(2);
%probfun = @(x_,y_) mvnpdf([x_(:), repmat(y_, size(x_(:)))], mu, sigma); % adapted to dblquadl
probfun = @(x_,y_) reshape(mvnpdf([x_(:), y_(:)], mu, sigma), size(x_)); % adapted to integral2

% Circular domain
rad = 0.5;
x_grid = -1:0.01:1;
y_grid = -1:0.01:1;
x_int = -rad:0.01:rad;

borderfun_pos = @(x_)real( sqrt(rad^2 - x_.^2));
borderfun_neg = @(x_)real(-sqrt(rad^2 - x_.^2));

[X,Y] = meshgrid(x_grid, y_grid);



prob = mvnpdf([X(:), Y(:)], mu, sigma);
prob = reshape(prob, length(x_grid), length(y_grid));
contour(X, Y, prob); hold on;
xlim([min(x_grid), max(x_grid)]);

phi = linspace(0,2*pi,1000); phi = phi(1:end-1);
x_circle = rad * cos(phi); y_circle = rad * sin(phi);
plot(x_int, borderfun_pos(x_int), 'b--');
plot(x_int, borderfun_neg(x_int), 'b--');
axis equal


% s = warning('off','MATLAB:quadl:MinStepSize');
% tol = [];
% % Inner integral (positive semicircle:
% inner_fun_pos = @(x_) quadl(@(y_) probfun(x_,y_), 0, borderfun_pos(x_));
% inner_fun_neg = @(x_) quadl(@(y_) probfun(x_,y_), 0, borderfun_pos(x_));

res = integral2(probfun, min(x_int), max(x_int), borderfun_neg, borderfun_pos);

