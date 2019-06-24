addpath('toolbox_signal/');
addpath('toolbox_general/');
addpath('toolbox_graph/');

n = 256; % size of the image
sigma = n/8; % width of the bumps
[Y,X] = meshgrid(1:n,1:n);
x = n/4; y = n/4;
M = exp( -( (X-x).^2 + (Y-y).^2 )/(2*sigma^2) );
x = 3*n/4; y = 3*n/4;
M = M + exp( -( (X-x).^2 + (Y-y).^2 )/(2*sigma^2) );

%W = rescale(M,3*1e-1,1);
W = ones(n);
%W(sqrt((X-n/2).^2 + (Y-n/2).^2) > 124) = 0;

figure; 
imageplot(M, [], 1,2,1);
imageplot(1./W, [], 1,2,2);

%% Farthest point sampling
profile on;
np = 1000;
tic; times = zeros(np,1);
vertices = zeros(2,np); vertices(:,1) = [89;37];
for i = 2:np
    D = perform_fast_marching(1./W, vertices(:, 1:i-1));
    [tmp,ind] = max(D(:));
    [x,y] = ind2sub([n n],ind);
    vertices(:,i) = [x;y];
    times(i) = toc;
end
%vertices = fp_sampling_wrapper( np, W, [1;1] );
profile report

%% Plotting
figure;
subplot(1,2,1);
hold on;
imageplot(W, 'Metric'); axis ij;
for i = 1:np
    %plot(vertices(2,:), vertices(1,:), 'r*');
    text(vertices(2,i), vertices(1,i), num2str(i));
end
subplot(1,2,2);
hold on;
imageplot( perform_hist_eq(D, 'linear'), 'Distance'); axis ij;
for i = 1:np
    %plot(vertices(2,:), vertices(1,:), 'r*');
    text(vertices(2,i), vertices(1,i), num2str(i));
end
colormap jet(256);

%%
figure; 
plot(sqrt(sum((vertices - [vertices(:, 2:end), [nan;nan]]).^2)));
%distmat = pdist2(vertices', vertices');
%distmat(distmat == 0) = inf;
%plot(min(distmat));


