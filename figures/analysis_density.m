element_a = 1; %mm^2
element_r = sqrt(element_a/pi);
valid_r = Rmuscle - element_r;

xg = linspace(-Rmuscle, Rmuscle, 100);
yg = linspace(-Rmuscle, Rmuscle, 100);
[Xg, Yg] = meshgrid(xg, yg);
valid_map = double(sqrt(Xg.^2 + Yg.^2) <= valid_r);
valid_map(valid_map == 0) = nan; 

%% Muscle fiber density in an element area
mf_density_result = zeros(size(Xg));
for i = 1:numel(xg)
    for j = 1:numel(yg)
        if sqrt(Xg(i,j).^2 + Yg(i,j).^2) < Rmuscle
            in_element = sqrt((mu_pool.mf_centers(:,1) - Xg(i,j)).^2 + (mu_pool.mf_centers(:,2) - Yg(i,j)).^2) < element_r;
            mf_density_result(i,j) = sum(in_element)/element_a;
        else
            mf_density_result(i,j) = nan;
        end
    end
end

figure; surf(Xg, Yg, mf_density_result .* valid_map); hold on;
%plot3(muscle_border(:,1), muscle_border(:,2), zeros(size(muscle_border(:,1))), 'r-');
%plot3(muscle_border(:,1)*(Rmuscle - element_r)/Rmuscle, muscle_border(:,2)*(Rmuscle - element_r)/Rmuscle, zeros(size(muscle_border(:,1))), 'k--');
%title('MF density over the muscle cross-section');
zlabel('Density of muscle fibers');
colorbar;
zlim([0, inf]);
figure2page;

%% Motor unit representation in an area element
mu_density_result = nan(size(Xg));

for i = 1:numel(xg)
    for j = 1:numel(yg)
        in_element = sqrt((mu_pool.mf_centers(:,1) - Xg(i,j)).^2 + (mu_pool.mf_centers(:,2) - Yg(i,j)).^2) < element_r;
        in_element = mu_pool.assignment(in_element);
        mu_density_result(i,j) = numel(unique(in_element));
    end
end

figure; 
mu_density_result(mu_density_result==0) = nan;
histogram(mu_density_result(:));
%title('Histogram of MU density over the muscle cross-section');

%%
mu_density_result = imgaussfilt(mu_density_result, 1);
figure; surf(Xg, Yg, mu_density_result .* valid_map); hold on;
%plot3(muscle_border(:,1)*(Rmuscle - element_r)/Rmuscle, muscle_border(:,2)*(Rmuscle - element_r)/Rmuscle, zeros(size(muscle_border(:,1))), 'k--');
zlim([0, inf]);
zlabel('Number of MNs present in a mm$^2$');
%title('MU density over the muscle cross-section');
colorbar;

clear mu_density_result mf_density_result


%% SD of fiber diameter in an element area
sd_diameter_result = zeros(size(Xg));
for i = 1:numel(xg)
    for j = 1:numel(yg)
        if sqrt(Xg(i,j).^2 + Yg(i,j).^2) < Rmuscle
            in_element = sqrt((mu_pool.mf_centers(:,1) - Xg(i,j)).^2 + (mu_pool.mf_centers(:,2) - Yg(i,j)).^2) < element_r;
            sd_diameter_result(i,j) = mean(mu_pool.mf_diameters(in_element)) * 1e6;%sum(in_element)/element_a;
        else
            sd_diameter_result(i,j) = nan;
        end
    end
end

figure; 
%surf(Xg, Yg, sd_diameter_result .* valid_map); hold on;
surf(Xg, Yg, sd_diameter_result .* valid_map); hold on;
%plot3(muscle_border(:,1), muscle_border(:,2), zeros(size(muscle_border(:,1))), 'r-');
%plot3(muscle_border(:,1)*(Rmuscle - element_r)/Rmuscle, muscle_border(:,2)*(Rmuscle - element_r)/Rmuscle, zeros(size(muscle_border(:,1))), 'k--');
%title('MF density over the muscle cross-section');
%zlabel('Mean fibers'' diameters in 1 $mm^2$ area, $\mu m$');
%title('Mean fibers'' diameters in 1 $mm^2$ area, $\mu m$');
xlabel('Muscle X axis, mm'); ylabel('Muscle Y axis, mm');
colorbar;
zlim([0, inf]);
figure2page;

clear valid_r valid_map xg yg Xg Yg in_element i j element_r element_a