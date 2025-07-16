function [bf_matrix, angle_vec] = fractionalBeamScan(rx_data, fs, fc, interpFactor, angle_range)
% fractionalBeamScan
% 对多个方向进行分数延迟波束形成
%
% 输入:
%   rx_data      : N×M 数据矩阵（时间 × 阵元）
%   fs           : 原始采样率 (Hz)
%   fc           : 中心频率 (Hz)
%   interpFactor : 插值倍数
%   angle_range  : 扫描角度向量 (单位：度)
%
% 输出:
%   bf_matrix    : L×K 矩阵，L为插值后采样点数，K为角度数
%   angle_vec    : 与列对应的扫描角度

    % 参数初始化
    c = 1500;
    [nSamples, Nelem] = size(rx_data);
    lambda = c / fc;
    d = lambda / 2;
    elem_pos = ((0:Nelem-1) - (Nelem-1)/2) * d;

    % 插值
    fs_interp = fs * interpFactor;
    rx_data_interp = resample(rx_data, fs_interp, fs);
    nInterp = size(rx_data_interp, 1);

    % 输出初始化
    angle_vec = angle_range(:)';
    nAngles = length(angle_vec);
    bf_matrix = zeros(nInterp, nAngles);

    % 遍历角度执行波束形成
    for k = 1:nAngles
        theta_rad = angle_vec(k) * pi / 180;
        steer_delays = elem_pos * sin(theta_rad) / c;
        steer_delays_interp = steer_delays * fs_interp;

        % 分数延迟求和
        bf = zeros(nInterp, 1);
        for m = 1:Nelem
            delay = steer_delays_interp(m);
            bf = bf + fracdelay(rx_data_interp(:, m), delay);
        end
        bf_matrix(:, k) = bf / Nelem;
    end
end

%% 内部函数：分数延迟
function y = fracdelay(x, d)
    N = 64;
    n = -N:N;
    h = sinc(n - d);
    h = h .* hamming(length(h))';
    h = h / sum(h);
    y = conv(x, h, 'same');
end
