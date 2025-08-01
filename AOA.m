%% aoa_ekf_multitrack_final.m

% Final Balanced EKF: High availability + low error
clear; close all; clc;
load('LNSM_Project_Data.mat');

%% Configuration

params.dt       = 0.1;
params.A        = [1 0 params.dt 0; 0 1 0 params.dt; 0 0 1 0; 0 0 0 1];
params.Q        = diag([0.01, 0.01, 0.1, 0.1]);  % smooth dynamics
params.sigma0   = 0.2;
params.R0       = params.sigma0^2;
params.P0       = eye(4) * 10;
params.MIN_UPDATES = 5;
params.SCALE_CLIP = [0.8, 2.0];  % moderate auto-tuning

% Constants
AP_pos = AP(1:2,:);
yaw    = APyaw;

track_names = { 'Track 1: with obstacle', 'Track 2: straight, no obstacle', 'Track 3: straight, noisy' };

%% Run for each track

for track = 1:3
    fprintf('\n===== %s =====\n', track_names{track});
    aoa_raw = AoA{track}(1:10,:);
    gt2d = ground_truth{track}(1:2,:);
    azg = wrapToPi(aoa_raw + yaw');  % no hard filtering

    N = size(azg,2);
    t = (0:N-1)*params.dt;

    [~, ~, nis1, ~] = run_ekf(params, AP_pos, azg);
    R_tuned = tune_R(nis1, azg, params);
    [x_est, ~, nis2, valid] = run_ekf(params, AP_pos, azg, R_tuned);

    [metrics, errs] = compute_metrics(x_est, gt2d, azg, nis2, valid, t);
    print_metrics(metrics);
    plot_results(x_est, gt2d, t, errs, nis2, valid, azg, track_names{track});
end

%% FUNCTIONS

function R = tune_R(nis, azg, p)
    upd = ~isnan(nis);
    if sum(upd) < p.MIN_UPDATES
        warning('Not enough updates for reliable tuning.');
        R = p.R0;
        return;
    end
    m_counts = sum(~isnan(azg), 1);
    scale = mean(nis(upd) ./ m_counts(upd));
    scale = min(max(scale, p.SCALE_CLIP(1)), p.SCALE_CLIP(2));
    R = p.R0 * scale;
    fprintf('Auto-tuned sigma^2 = %.4f (scale = %.2f)\n', R, scale);
end

function [metrics, errs] = compute_metrics(x, gt, azg, nis, valid, t)
    est = x(1:2, valid);
    gt_valid = gt(:, valid);
    errs = sqrt(sum((est - gt_valid).^2));
    metrics.MAE = mean(errs);
    metrics.RMSE = sqrt(mean(errs.^2));
    metrics.Q95 = prctile(errs, 95);
    metrics.Availability = 100 * sum(valid) / length(valid);
    m_counts = sum(~isnan(azg),1);
    chi2th = chi2inv(0.95, m_counts);
    goodNIS = valid & (nis <= chi2th);
    metrics.Reliability = 100 * sum(goodNIS) / sum(valid);
end

function print_metrics(m)
    fprintf('MAE = %.3f m\n', m.MAE);
    fprintf('RMSE = %.3f m\n', m.RMSE);
    fprintf('95%% quantile = %.3f m\n', m.Q95);
    fprintf('Availability = %.1f %%\n', m.Availability);
    fprintf('Reliability = %.1f %%\n', m.Reliability);
end

function plot_results(x, gt, t, errs, nis, valid, azg, name)
    m_counts = sum(~isnan(azg),1);
    chi2th = chi2inv(0.95, m_counts);

    figure('Name', name);
    sgtitle(name);

    subplot(2,2,1); hold on; grid on; axis equal;
    plot(gt(1,:), gt(2,:), 'b-', 'LineWidth', 1.4);
    plot(x(1,:), x(2,:), 'r-', 'LineWidth', 1.4);
    title('Trajectory'); xlabel('X [m]'); ylabel('Y [m]'); legend('GPS','EKF');

    subplot(2,2,2);
    plot(t(valid), errs, 'k', 'LineWidth', 1.2); grid on;
    title('Error vs Time'); xlabel('Time [s]'); ylabel('Error [m]');

    subplot(2,2,3);
    [f,x_e] = ecdf(errs);
    plot(x_e, f, 'm', 'LineWidth', 1.4); grid on;
    title('Error CDF'); xlabel('Error [m]'); ylabel('CDF');

    subplot(2,2,4);
    plot(t(valid), nis(valid), 'g', 'LineWidth', 1.2); hold on;
    plot(t(valid), chi2th(valid), '--k'); grid on;
    title('NIS vs Time'); xlabel('Time [s]'); ylabel('NIS');
    legend('NIS','\chi^2_{0.95}','Interpreter','tex');

    vx = x(3,valid); vy = x(4,valid);
    speed = sqrt(vx.^2 + vy.^2);
    figure('Name', [name ' Speed']);
    plot(t(valid), speed, 'b-', 'LineWidth', 1.4); grid on;
    title('Estimated Speed'); xlabel('Time [s]'); ylabel('Speed [m/s]');
end

function [x_est,P,nis_list,valid_mask] = run_ekf(p, AP, azg, R_scalar)
    if nargin < 4, R_scalar = p.R0; end
    N = size(azg,2); x = NaN(4,1); P = p.P0;
    x_est = NaN(4,N); nis_list = nan(1,N); valid_mask = false(1,N);

    for t = 1:N
        if t > 1, x = p.A*x; P = p.A*P*p.A' + p.Q; end

        if isnan(x(1))
            p0 = aoa_ls_initial(AP, azg(:,t));
            if any(isnan(p0)), continue; end
            x(1:2) = p0; x(3:4) = 0; x_est(:,t)=x; valid_mask(t)=true; continue;
        end

        idx = find(~isnan(azg(:,t))); m = numel(idx);
        if m < 2, continue; end

        z = azg(idx,t); R = R_scalar * eye(m); h = zeros(m,1); Hm = zeros(m,4);
        for k = 1:m
            i = idx(k); dx = x(1)-AP(1,i); dy = x(2)-AP(2,i); d2 = dx^2+dy^2;
            h(k) = atan2(dy,dx); Hm(k,1) = -dy/d2; Hm(k,2) = dx/d2;
        end

        nu = wrapToPi(z - h); S = Hm*P*Hm' + R; K = P*Hm'/S;
        x = x + K*nu; P = (eye(4)-K*Hm)*P;
        x_est(:,t) = x; nis_list(t) = nu'/S*nu; valid_mask(t)=true;
    end
end

function pos = aoa_ls_initial(AP, az)
    idx = find(~isnan(az)); if numel(idx)<3, pos=[NaN;NaN]; return; end
    a = wrapToPi(az(idx)); med = median(a);
    dev = abs(wrapToPi(a - med)); w = exp(-4*dev);
    A = zeros(numel(idx),2); b = zeros(numel(idx),1);
    for k = 1:numel(idx)
        i = idx(k); dir = [cos(az(i)); sin(az(i))]; norm_v = [-dir(2); dir(1)];
        A(k,:) = w(k)*norm_v'; b(k) = w(k)*(norm_v'*AP(:,i));
    end
    if rcond(A'*A)<1e-6, pos=[NaN;NaN]; else pos = (A'*A) \ (A'*b); end
end