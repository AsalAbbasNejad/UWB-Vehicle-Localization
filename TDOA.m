clc; clear; close all;
load('LNSM_Project_Data.mat');

Ts = 0.1;
chi2_th = chi2inv(0.95, 2); % Reliability threshold for 2 DoF

%% ===== Track 1: with obstacle =====
tdoa = TDoA{1};
gt = ground_truth{1};
AP = AP;

positions = NaN(2, size(tdoa, 2));
NIS_list = NaN(1, size(tdoa, 2));

for t = 1:size(tdoa, 2)
    master = tdoa(10, t);
    if isnan(master), continue; end

    tdoa_col = tdoa(1:9, t);
    tdoa_idx = find(~isnan(tdoa_col));
    if numel(tdoa_idx) < 3, continue; end

    Ai = AP(1:2, tdoa_idx);
    Am = AP(1:2, master);
    rho = tdoa_col(tdoa_idx);

    fun = @(p) vecnorm(p - Ai, 2, 1) - vecnorm(p - Am, 2, 1) - rho';
    pos0 = mean(Ai, 2);

    [pos, resnorm] = lsqnonlin(fun, pos0, [], [], ...
        optimoptions('lsqnonlin','Display','off','MaxIter',50,'TolFun',1e-4));

    if resnorm < 20
        positions(:, t) = pos;
    end
end

for i = 1:2
    positions(i,:) = hampel(positions(i,:), 5, 3);
end

A_kf = [1 0 Ts 0; 0 1 0 Ts; 0 0 1 0; 0 0 0 1];
H_kf = [1 0 0 0; 0 1 0 0];
Q_kf = 0.05 * eye(4);
R_kf = 2.0 * eye(2);
P_kf = 10 * eye(4);
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

valid = all(~isnan(positions_kf), 1) & all(~isnan(gt(1:2,:)), 1);
est_valid = positions_kf(:, valid);
gt_valid = gt(1:2, valid);
errors = sqrt(sum((est_valid - gt_valid).^2, 1));

MAE = mean(errors);
RMSE = sqrt(mean(errors.^2));
Q95 = prctile(errors, 95);
availability = sum(valid) / size(positions_kf, 2) * 100;
reliability = 100 * sum(NIS_list(valid) <= chi2_th) / sum(valid);

fprintf('\n===== Track 1: with obstacle =====\n');
fprintf('MAE: %.2f m\n', MAE);
fprintf('RMSE: %.2f m\n', RMSE);
fprintf('Q95: %.2f m\n', Q95);
fprintf('Availability: %.2f %%\n', availability);
fprintf('Reliability: %.2f %%\n', reliability);

t_sec = (0:size(positions_kf,2)-1) * Ts;
figure('Name','Track 1: with obstacle','NumberTitle','off');
sgtitle('TDOA-only Improved — Track 1');

subplot(2,2,1);
plot(gt(1,valid), gt(2,valid), 'b', 'LineWidth',1.5); hold on;
plot(est_valid(1,:), est_valid(2,:), 'r');
legend('Ground Truth','TDOA Estimate','Location','best');
xlabel('X [m]'); ylabel('Y [m]'); title('Trajectory');
axis equal; grid on;

subplot(2,2,2);
plot(t_sec(valid), errors, 'k', 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Error [m]');
title('Position Error vs Time');

subplot(2,2,3);
[f_e, x_e] = ecdf(errors);
plot(x_e, f_e, 'm', 'LineWidth', 1.4); grid on;
xlabel('Error [m]'); ylabel('CDF');
title('Error CDF');

subplot(2,2,4);
vel_est = diff(est_valid,1,2) / Ts;
vel_mag = sqrt(sum(vel_est.^2, 1));
t_v = t_sec(2:end);
min_len = min(length(t_v), length(vel_mag));
plot(t_v(1:min_len), vel_mag(1:min_len), 'b', 'LineWidth', 1.4); grid on;
xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Estimated Speed');

%% ===== Track 2 and 3 Combined =====

track_names = { 'Track 2: straight', 'Track 3: noisy' };
for track = 2:3
    fprintf('\n===== %s =====\n', track_names{track - 1});

    tdoa = TDoA{track};
    gt = ground_truth{track};
    AP = AP;

    positions = NaN(2, size(tdoa, 2));
    NIS_list = NaN(1, size(tdoa, 2));

    for t = 1:size(tdoa, 2)
        master = tdoa(10, t);
        if isnan(master), continue; end

        tdoa_col = tdoa(1:9, t);
        tdoa_idx = find(~isnan(tdoa_col));
        tdoa_idx(tdoa_idx == master) = [];
        if numel(tdoa_idx) < 3, continue; end

        Ai = AP(1:2, tdoa_idx);
        Am = AP(1:2, master);
        rho = tdoa_col(tdoa_idx);

        fun = @(p) vecnorm(p - Ai, 2, 1) - vecnorm(p - Am, 2, 1) - rho';
        pos0 = mean(Ai, 2);

        [pos, resnorm] = lsqnonlin(fun, pos0, [], [], ...
            optimoptions('lsqnonlin','Display','off'));

        if resnorm < 10
            positions(:, t) = pos;
        end
    end

    for i = 1:2
        positions(i,:) = hampel(positions(i,:), 7, 2);
    end

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

    valid = all(~isnan(positions_kf), 1) & all(~isnan(gt(1:2,:)), 1);
    est_valid = positions_kf(:, valid);
    gt_valid = gt(1:2, valid);
    errors = sqrt(sum((est_valid - gt_valid).^2, 1));

    MAE = mean(errors);
    RMSE = sqrt(mean(errors.^2));
    Q95 = prctile(errors, 95);
    availability = sum(valid) / size(positions_kf, 2) * 100;
    reliability = 100 * sum(NIS_list(valid) <= chi2_th) / sum(valid);

    fprintf('MAE: %.2f m\n', MAE);
    fprintf('RMSE: %.2f m\n', RMSE);
    fprintf('Q95: %.2f m\n', Q95);
    fprintf('Availability: %.2f %%\n', availability);
    fprintf('Reliability: %.2f %%\n', reliability);

    t_sec = (0:size(positions_kf, 2)-1) * Ts;
    figure('Name',track_names{track - 1},'NumberTitle','off');
    sgtitle(['TDOA-only — ' track_names{track - 1}]);

    subplot(2,2,1);
    plot(gt(1,valid), gt(2,valid), 'b', 'LineWidth',1.5); hold on;
    plot(est_valid(1,:), est_valid(2,:), 'r');
    legend('Ground Truth','TDOA Estimate','Location','best');
    xlabel('X [m]'); ylabel('Y [m]'); title('Trajectory');
    axis equal; grid on;

    subplot(2,2,2);
    plot(t_sec(valid), errors, 'k', 'LineWidth', 1.2); grid on;
    xlabel('Time [s]'); ylabel('Error [m]');
    title('Position Error vs Time');

    subplot(2,2,3);
    [f_e, x_e] = ecdf(errors);
    plot(x_e, f_e, 'm', 'LineWidth', 1.4); grid on;
    xlabel('Error [m]'); ylabel('CDF');
    title('Error CDF');

    subplot(2,2,4);
    vel_est = diff(est_valid,1,2) / Ts;
    vel_mag = sqrt(sum(vel_est.^2, 1));
    t_v = t_sec(2:end);
    min_len = min(length(t_v), length(vel_mag));
    plot(t_v(1:min_len), vel_mag(1:min_len), 'b', 'LineWidth', 1.4); grid on;
    xlabel('Time [s]'); ylabel('Speed [m/s]');
    title('Estimated Speed');
end
