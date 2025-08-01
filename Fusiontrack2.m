%% TDOA + AOA Final with Reliability

clc; clear; close all;
load('LNSM_Project_Data.mat');

track = 2;
tdoa = TDoA{track};
aoa = AoA{track};
yaw = APyaw;
gt = ground_truth{track};
AP = AP;
Ts = 0.1;

az_global = aoa(1:10,:) + yaw';
positions = NaN(2, size(tdoa, 2));
NIS_list = NaN(1, size(tdoa,2));  % For reliability

for t = 1:size(tdoa, 2)
    master = tdoa(10, t);
    if isnan(master), continue; end

    tdoa_col = tdoa(1:9, t);
    az_col = az_global(:, t);

    tdoa_idx = find(~isnan(tdoa_col));
    aoa_idx = find(~isnan(az_col));
    common_idx = intersect(tdoa_idx, aoa_idx);

    % --- Primary Fusion: if enough TDOA + AOA available ---
    if numel(common_idx) >= 2 && numel(aoa_idx) >= 1
        Ai = AP(1:2, common_idx);
        Am = AP(1:2, master);
        rho = tdoa_col(common_idx);
        az = az_col(common_idx);

        fun = @(p) [
            vecnorm(p - Ai, 2, 1) - vecnorm(p - Am, 2, 1) - rho';
            wrapToPi(atan2(p(2) - Ai(2,:), p(1) - Ai(1,:)) - az')
        ];

        pos0 = mean(Ai, 2);
        [pos, resnorm] = lsqnonlin(fun, pos0, [], [], optimoptions('lsqnonlin','Display','off'));

        if resnorm < 9
            positions(:, t) = pos;
            NIS_list(t) = resnorm;  % Save squared residual norm as NIS
            continue;
        end
    end

    % --- Fallback to AOA-only ---
    if numel(aoa_idx) >= 2
        A_aoa = AP(1:2, aoa_idx);
        az = az_col(aoa_idx);

        fun_aoa = @(p) wrapToPi(atan2(p(2) - A_aoa(2,:), p(1) - A_aoa(1,:)) - az');
        pos0 = mean(A_aoa, 2);

        [pos_aoa, res_aoa] = lsqnonlin(fun_aoa, pos0, [], [], optimoptions('lsqnonlin','Display','off'));

        if res_aoa < 5
            positions(:, t) = pos_aoa;
            NIS_list(t) = res_aoa;  % Store AOA-only residual norm
        end
    end
end

valid = all(~isnan(positions), 1) & all(~isnan(gt(1:2,:)), 1);
est_valid = positions(:, valid);
gt_valid = gt(1:2, valid);

errors = sqrt(sum((est_valid - gt_valid).^2, 1));
MAE = mean(errors);
RMSE = sqrt(mean(errors.^2));
availability = sum(valid) / size(positions, 2) * 100;

% Reliability threshold: DoF = number of constraints = 2D (x,y)
chi2_th = chi2inv(0.95, 2);
reliability = 100 * sum(NIS_list(valid) <= chi2_th) / sum(valid);

fprintf('--- Final TDOA + AOA Fusion ---\n');
fprintf('MAE: %.2f m\n', MAE);
fprintf('RMSE: %.2f m\n', RMSE);
fprintf('Availability: %.2f %%\n', availability);
fprintf('Reliability: %.2f %%\n', reliability);

t_sec = (0:size(tdoa,2)-1) * Ts;

figure;
plot(gt(1,valid), gt(2,valid), 'b', 'LineWidth', 1.5); hold on;
plot(est_valid(1,:), est_valid(2,:), 'r');
legend('Ground Truth', 'TDOA+AOA Estimate');
xlabel('X [m]'); ylabel('Y [m]');
title('2D Localization using TDOA + AOA');
axis equal; grid on;

figure;
plot(t_sec(valid), errors, 'g', 'LineWidth', 1.4); grid on;
xlabel('Time [s]');
ylabel('Localization Error [m]');
title('TDOA + AOA - Error over Time');

vel_est = diff(est_valid,1,2) / Ts;
vel_mag = sqrt(sum(vel_est.^2, 1));
t_v = t_sec(valid);
t_v = t_v(2:end);  % match velocity vector length

figure;
plot(t_v, vel_mag, 'LineWidth', 1.5); grid on;
xlabel('Time [s]');
ylabel('Estimated Velocity [m/s]');
title('Estimated Speed over Time - [TDOA + AOA]');

figure;
plot(t_sec(valid), NIS_list(valid), 'm', 'LineWidth', 1.4); hold on;
yline(chi2_th, '--r', '\chi^2_{0.95}');
xlabel('Time [s]');
ylabel('NIS Value');
title('NIS vs Time — TDOA + AOA');
grid on;
