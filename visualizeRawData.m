% visualizeRawDataStandalone.m
% Standalone script for Task 1 & 3: Load, visualize raw UWB (TDoA) and AoA data,
%  static 3D rays, and animate 3D AoA rays over time

clear; close all; clc;

%% Settings
dataFile     = 'LNSM_Project_Data.mat';   % .mat file containing AP, AoA, TDoA, GPS
datasetIndex = 1;                          % select dataset 1, 2, or 3
Fs           = 10;                         % sample rate [Hz]
rayLen       = 50;                         % length of each AoA ray [m]
tPlotStatic  = 1;                          % time index for static 3D view (1<=tPlotStatic<=N)

%% Load data
data = loadData(dataFile, datasetIndex);
N    = size(data.rawTDoA, 2);
time = (0:N-1) / Fs;
APpos = data.APpos;                        % 10×3
rawAoA = data.rawAoA;                      % 20×N
GT     = data.GPS;                         % 3×N

%% Raw TDoA channels (1–9)

figure('Name','Raw TDoA Channels','NumberTitle','off'); % First Track
tdoaRaw = data.rawTDoA(1:9, :);            % rows 1–9
plot(time, tdoaRaw', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('TDoA [m]');
title('Raw TDoA Measurements (channels 1–9)'); grid on;

%% Raw AoA azimuths (rows 1–10)
figure('Name','Raw AoA Azimuths','NumberTitle','off');
aoaAz = rawAoA(1:10, :);
plot(time, rad2deg(aoaAz)', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('Azimuth [°]');
title('Raw AoA Azimuth Measurements (per AP)');
legend(arrayfun(@(i) sprintf('AP%d',i), 1:10, 'UniformOutput', false), 'Location', 'best'); grid on;

%% Raw AoA elevations (rows 11–20)
figure('Name','Raw AoA Elevations','NumberTitle','off');
aoaEl = rawAoA(11:20, :);
plot(time, rad2deg(aoaEl)', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('Elevation [°]');
title('Raw AoA Elevation Measurements (per AP)');
legend(arrayfun(@(i) sprintf('AP%d',i), 1:10, 'UniformOutput', false), 'Location', 'best'); grid on;

%% Static 3D AoA Rays Visualization (Task 3)

tPlot = tPlotStatic;
az = rawAoA(1:10, tPlot);
el = rawAoA(11:20, tPlot);
dirsStatic = [cos(el).*cos(az), cos(el).*sin(az), sin(el)];

figure('Name','3D AoA Rays (Static)','NumberTitle','off');
hold on; grid on; axis equal;
scatter3(APpos(:,1), APpos(:,2), APpos(:,3), 100, 'k', 'filled');
text(APpos(:,1), APpos(:,2), APpos(:,3), arrayfun(@(i)sprintf('AP%d',i),1:10,'Uni',0), ...
     'VerticalAlignment','bottom','HorizontalAlignment','right');
for i = 1:10
    P0 = APpos(i,:);
    P1 = P0 + rayLen * dirsStatic(i,:);
    plot3([P0(1),P1(1)], [P0(2),P1(2)], [P0(3),P1(3)], 'r-', 'LineWidth',1.5);
end
plot3(GT(1,:), GT(2,:), GT(3,:), 'b-', 'LineWidth',1.5);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(sprintf('3D AoA Rays at t = %d (%.1f s)', tPlot, (tPlot-1)/Fs));
legend('AP positions','AoA rays','GPS RTK GT','Location','best');
view(3); rotate3d on;

%% Animate 3D AoA Rays Over Time
figure('Name','3D AoA Rays Animation','NumberTitle','off');
hold on; grid on; axis equal;
% Plot AP positions once
t = scatter3(APpos(:,1), APpos(:,2), APpos(:,3), 100, 'k', 'filled');
text(APpos(:,1), APpos(:,2), APpos(:,3), arrayfun(@(i)sprintf('AP%d',i),1:10,'Uni',0), ...
     'VerticalAlignment','bottom','HorizontalAlignment','right');
% Initialize ray and GT handles
rayH = gobjects(10,1);
for i = 1:10
    rayH(i) = plot3(NaN,NaN,NaN,'r-','LineWidth',1);
end
gtH = plot3(NaN,NaN,NaN,'b-','LineWidth',1.5);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
view(3); rotate3d on;

title('3D AoA Rays Animation');
for k = 1:N
    azk = rawAoA(1:10, k);
    elk = rawAoA(11:20, k);
    dirs = [cos(elk).*cos(azk), cos(elk).*sin(azk), sin(elk)];
    % update each ray
    for i = 1:10
        P0 = APpos(i,:);
        P1 = P0 + rayLen * dirs(i,:);
        set(rayH(i), 'XData', [P0(1),P1(1)], ...
                     'YData', [P0(2),P1(2)], ...
                     'ZData', [P0(3),P1(3)]);
    end
    % update GT trace
    set(gtH, 'XData', GT(1,1:k), 'YData', GT(2,1:k), 'ZData', GT(3,1:k));
    title(sprintf('Frame %d/%d, Time = %.1f s', k, N, (k-1)/Fs));
    drawnow;
end

%% Local function: loadData
function data = loadData(filename, idx)
    tmp          = load(filename);
    data.APpos   = tmp.AP';               % 10×3
    data.APyaw   = tmp.APyaw';            % 10×1 (not used here)
    data.APpitch = tmp.APpitch';          % 10×1 (not used here)
    data.rawTDoA = tmp.TDoA{idx};         % 10×N
    data.rawAoA  = tmp.AoA{idx};          % 20×N
    data.GPS     = tmp.ground_truth{idx}; % 3×N
end
