function sig = createLFM(BW, T, fs, f0, sweepDir)
% createLFM  生成实信号 LFM 脉冲（幅值归一化）
%
% sig = createLFM(BW, T, fs, f0, sweepDir)
%
% BW        : 调频带宽 (Hz)，如 1e5
% T         : 脉宽 (s)，   如 5e-3
% fs        : 采样率 (Hz)，如 2e6
% f0        : 起始频率 (Hz)，可省略；默认 sweep 以 (f0, f0+BW)
% sweepDir  : 'up' 或 'down'，可省略；默认 'up'
%
% 返回 sig   : 1×N 实向量，幅值已归一化到 ±1

    arguments
        BW   (1,1) double {mustBePositive}
        T    (1,1) double {mustBePositive}
        fs   (1,1) double {mustBePositive}
        f0   (1,1) double {mustBeNonnegative} = 0
        sweepDir (1,1) string {mustBeMember(sweepDir,["up","down"])} = "up"
    end

    t  = (0:1/fs:T-1/fs);          % 时间轴
    k  = BW/T;                     % 调频斜率
    if sweepDir == "down", k = -k; end
    instPhase = 2*pi*(f0.*t + 0.5*k.*t.^2);
    sig = cos(instPhase);

    % 幅值归一化到 ±1
    sig = sig ./ max(abs(sig));
end
