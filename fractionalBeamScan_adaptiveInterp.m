function [bf_matrix, angle_vec] = fractionalBeamScan_adaptiveInterp(rx_data, fs, fc, angle_range, interp_eps)
% 自适应插值倍数的分数延迟波束形成
%
% 输入:
%   rx_data     : N×M 输入数据（时间 × 阵元）
%   fs          : 采样率 (Hz)
%   fc          : 中心频率 (Hz)
%   angle_range : 扫描角度向量（单位：度）
%   interp_eps  : 允许的最大分数延迟误差（单位：采样点，推荐0.1）
%
% 输出:
%   bf_matrix   : L×K 矩阵，L为插值后采样点数，K为方向数
%   angle_vec   : K×1 扫描角度（deg）

    if nargin < 5
        interp_eps = 0.1; % 默认允许分数延迟误差为0.1采样点
    end

    % 参数初始化
    c = 1500;
    [nSamples, Nelem] = size(rx_data);
    lambda = c / fc;
    d = lambda / 2;
    elem_pos = ((0:Nelem-1) - (Nelem-1)/2) * d;

    angle_vec = angle_range(:)';
    nAngles = length(angle_vec);
    bf_matrix = cell(1, nAngles);

    % 针对每个角度，自适应计算最小插值倍数
    for k = 1:nAngles
        theta_rad = angle_vec(k) * pi / 180;
        delays = elem_pos * sin(theta_rad) / c;        % 秒
        delay_in_samples = delays * fs;

        frac_part = abs(delay_in_samples - round(delay_in_samples));
        max_frac = max(frac_part);

        % 自适应插值倍数：最小使误差小于 interp_eps
        interpFactor = ceil(1 / interp_eps * max_frac);
        if interpFactor < 1
            interpFactor = 1;
        elseif mod(interpFactor,2) == 1
            interpFactor = interpFactor + 1; % 保证偶数插值因子（FIR对称性）
        end

        % 插值数据
        fs_interp = fs * interpFactor;
        rx_interp = resample(rx_data, fs_interp, fs);  % [N_interp × Nelem]
        nInterp = size(rx_interp,1);

        steer_delay_interp = delays * fs_interp;

        % 波束形成
        bf = zeros(nInterp, 1);
        for m = 1:Nelem
            bf = bf + fracdelay(rx_interp(:, m), steer_delay_interp(m));
        end
        bf_matrix{k} = bf / Nelem;
    end

    % 转换 cell 为矩阵（按最小长度统一裁剪）
    minL = min(cellfun(@length, bf_matrix));
    bf_matrix = cell2mat(cellfun(@(x) x(1:minL), bf_matrix, 'UniformOutput', false));
end

%% 分数延迟函数（带窗 sinc 插值）
function y = fracdelay(x, d)
    N = 64;
    n = -N:N;
    h = sinc(n - d);
    h = h .* hamming(length(h))';
    h = h / sum(h);
    y = conv(x, h, 'same');
end
