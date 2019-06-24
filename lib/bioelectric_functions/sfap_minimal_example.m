L1 = 25; %mm
L2 = 25; %mm
z_endplate = 25; %mm

z_electrode = 40; %mm
r_electrode = 0.1; %mm

v = 4000; %mm/s

dz = 1e-2;
z = -10:dz:60;
dt = 1e-4;
t = 0:dt:10e-3;

h = get_elementary_current_response(z, z_electrode, r_electrode);
i = get_current_density(t, z, z_endplate, L1, L2, v);
response = h * i;

figure; plot(response); xlabel('z'); 
[T,Z] = meshgrid(t, z);
figure; mesh(Z,T,i); xlabel('z'); ylabel('Time');

