function bf_down = downsampleBeamData(bf_data, fs_orig, fs_target)
% 将波束输出从高采样率降采到目标采样率
%
% 输入:
%   bf_data   : N×K 波束数据（时间 × 角度）
%   fs_orig   : 原采样率（Hz，可能是插值后的采样率）
%   fs_target : 目标采样率（Hz，通常是原始 fs）
%
% 输出:
%   bf_down   : M×K 降采样后数据

    [N, K] = size(bf_data);
    [p, q] = rat(fs_target / fs_orig, 1e-6);   % 分数形式表示采样率比

    bf_down = zeros(ceil(N * p / q), K);
    for k = 1:K
        bf_down(:, k) = resample(bf_data(:, k), p, q);
    end
end
