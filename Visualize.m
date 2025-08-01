%% Display UWB Antenna Positions in 3D

figure;
scatter3(AP(1,:), AP(2,:), AP(3,:), 80, 'filled');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('UWB Antenna Positions');
grid on;

%% Ground Truth Trajectory (Free Track)

gt = ground_truth{2}; % tracks over time
figure;
plot3(gt(1,:), gt(2,:), gt(3,:), 'LineWidth', 2);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Ground Truth Trajectory (Free Track)');
grid on;

%% Global Azimuth - AP1 

aoa = AoA{1};             
az_raw = aoa(1:10, :);   
el_raw = aoa(11:20, :);
az_global = az_raw + APyaw';
el_global = el_raw + APpitch';                                                                                                                                                                                                              
plot(az_global(1, :)); title('Global Azimuth - AP1');

%% 

% Extract AOA
aoa = AoA{2};
az_raw = aoa(1:10,:);
el_raw = aoa(11:20,:);
az_global = az_raw + APyaw';
el_global = el_raw + APpitch';


% Plot Elevation of All Antennas
figure;
plot(rad2deg(el_global'));
xlabel('Sample Index');
ylabel('Elevation [deg]');
title('Global Elevation of All Antennas - Test 2');
grid on;
legend(arrayfun(@(i) sprintf('AP%d', i), 1:10, 'UniformOutput', false));


% Compare Azimuth and Elevation of AP1
figure;
subplot(2,1,1);
plot(rad2deg(az_global(1,:)));
xlabel('Sample Index');
ylabel('Azimuth [deg]');
title('Global Azimuth - AP1');
grid on;

subplot(2,1,2);
plot(rad2deg(el_global(1,:)));
xlabel('Sample Index');
ylabel('Elevation [deg]');
title('Global Elevation - AP1');
grid on;
