function bf_compressed = pulseCompress(bf_matrix, fs_interp, lfm_params)
% 对波束形成后的结果做脉冲压缩
%
% 输入:
%   bf_matrix   : N×K 波束输出（时间 × 角度）
%   fs_interp   : 插值后的采样率（Hz）
%   lfm_params  : struct 结构体，包含:
%                 .f0     : 起始频率 (Hz)
%                 .f1     : 结束频率 (Hz)
%                 .T_p    : 脉冲宽度 (s)
%
% 输出:
%   bf_compressed : N×K 脉压输出（时间 × 角度）

    % 构造匹配滤波器（LFM回波）
    t_p = 0:1/fs_interp:lfm_params.T_p - 1/fs_interp;
    lfm = chirp(t_p, lfm_params.f0, lfm_params.T_p, lfm_params.f1);
    lfm = lfm .* hamming(length(lfm))';  % 可选窗函数
    h = conj(flipud(lfm(:)));            % 匹配滤波器

    [nTime, nBeam] = size(bf_matrix);
    bf_compressed = zeros(nTime, nBeam);
    for k = 1:nBeam
        y = conv(bf_matrix(:, k), h, 'same');
        bf_compressed(:, k) = y;
    end
end
