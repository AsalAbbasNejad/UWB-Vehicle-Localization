clc; clear; close all;
load('LNSM_Project_Data.mat');

track = 1;
tdoa = TDoA{track};
aoa = AoA{track};
yaw = APyaw;
gt = ground_truth{track};
AP = AP;
Ts = 0.1;

az_global = wrapToPi(aoa(1:10,:) + yaw');
positions = NaN(2, size(tdoa, 2));
NIS_list = NaN(1, size(tdoa, 2));  % For reliability
chi2_th = chi2inv(0.95, 2);        % Threshold for 95% confidence (2 DoF)

for t = 1:size(tdoa, 2)
    master = tdoa(10, t);
    if isnan(master), continue; end

    tdoa_col = tdoa(1:9, t);
    az_col = az_global(:, t);
    tdoa_idx = find(~isnan(tdoa_col));
    aoa_idx = find(~isnan(az_col));
    common_idx = intersect(tdoa_idx, aoa_idx);

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
            continue;
        end
    end

    if numel(aoa_idx) >= 2
        A_aoa = AP(1:2, aoa_idx);
        az = az_col(aoa_idx);

        fun_aoa = @(p) wrapToPi(atan2(p(2) - A_aoa(2,:), p(1) - A_aoa(1,:)) - az');
        pos0 = mean(A_aoa, 2);
        [pos_aoa, res_aoa] = lsqnonlin(fun_aoa, pos0, [], [], optimoptions('lsqnonlin','Display','off'));
        if res_aoa < 5
            positions(:, t) = pos_aoa;
        end
    end
end

% Hampel Filter
for i = 1:2
    positions(i,:) = hampel(positions(i,:), 7, 2);
end

% Kalman Filter
A_kf = [1 0 Ts 0; 0 1 0 Ts; 0 0 1 0; 0 0 0 1];
H_kf = [1 0 0 0; 0 1 0 0];
Q_kf = 0.01 * eye(4);
R_kf = 1.0 * eye(2);
P_kf = eye(4);
x_kf = [NaN; NaN; 0; 0];

positions_kf = NaN(2, size(positions, 2));
for t = 1:size(positions, 2)
    meas = positions(:, t);
    if any(isnan(meas)), continue; end
    if isnan(x_kf(1))
        x_kf(1:2) = meas;
        positions_kf(:, t) = meas;
        continue;
    end

    x_pred = A_kf * x_kf;
    P_pred = A_kf * P_kf * A_kf' + Q_kf;
    y = meas - H_kf * x_pred;
    S = H_kf * P_pred * H_kf' + R_kf;
    K = P_pred * H_kf' / S;
    x_kf = x_pred + K * y;
    P_kf = (eye(4) - K * H_kf) * P_pred;

    positions_kf(:, t) = x_kf(1:2);
    NIS_list(t) = y' / S * y;
end

% Evaluation
valid = all(~isnan(positions_kf), 1) & all(~isnan(gt(1:2,:)), 1);
est_valid = positions_kf(:, valid);
gt_valid = gt(1:2, valid);

errors_fusion = sqrt(sum((est_valid - gt_valid).^2, 1));
MAE_Fusion = mean(errors_fusion);
RMSE_Fusion = sqrt(mean(errors_fusion.^2));
availability_Fusion = sum(valid) / size(positions, 2) * 100;
reliability_Fusion = 100 * sum(NIS_list(valid) <= chi2_th) / sum(valid);

% Output
fprintf('--- Final TDOA + AOA Fusion (Improved) ---\n');
fprintf('MAE: %.2f m\n', MAE_Fusion);
fprintf('RMSE: %.2f m\n', RMSE_Fusion);
fprintf('Availability: %.2f %%\n', availability_Fusion);
fprintf('Reliability: %.2f %%\n', reliability_Fusion);

% Plots
t_sec = (0:size(positions,2)-1) * Ts;
figure;
plot(gt(1,valid), gt(2,valid), 'b', 'LineWidth', 1.5); hold on;
plot(est_valid(1,:), est_valid(2,:), 'r');
legend('Ground Truth', 'TDOA+AOA Estimate');
xlabel('X [m]'); ylabel('Y [m]');
title('2D Localization using TDOA + AOA - Improved');
axis equal; grid on;

figure;
plot(t_sec(valid), errors_fusion, 'g', 'LineWidth', 1.4); grid on;
xlabel('Time [s]'); ylabel('Localization Error [m]');
title('TDOA + AOA - Error over Time');

vel_est = diff(est_valid,1,2) / Ts;
vel_mag = sqrt(sum(vel_est.^2, 1));
t_v = t_sec(2:end);
min_len = min(length(t_v), length(vel_mag));
figure;
plot(t_v(1:min_len), vel_mag(1:min_len), 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Estimated Velocity [m/s]');
title('Estimated Speed over Time - TDOA + AOA');
grid on;

% Optional: NIS Plot
figure;
valid_nis = valid & ~isnan(NIS_list);
plot(t_sec(valid_nis), NIS_list(valid_nis), 'm', 'LineWidth', 1.4); hold on;
yline(chi2_th, '--r', '\chi^2_{0.95}');
xlabel('Time [s]'); ylabel('NIS');
title('NIS vs Time — TDOA + AOA Fusion');
grid on;
