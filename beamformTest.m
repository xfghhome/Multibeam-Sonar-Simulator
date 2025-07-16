clear all;close all;clc
fs = 2e6;
fc = 900e3;
interpFactor = 8;
angle_range = -60:1:60;    % 扫描角度范围
eps = 0.1;                  % 插值误差容限（<1/10个采样点）
c_water = 1500;
lfm_params.f0 = 850e3;        % 起始频率
lfm_params.f1 = 950e3;        % 终止频率
lfm_params.T_p = 5e-3;        % 脉冲宽度 5 ms

for i = 1:11
    dataFileName = ['./OutputChannelData/ping_',num2str(i),'.dat'];
    pingData = loadPingData(dataFileName,'single');
    % 输入数据 rx_data: (50000×128)，采样率2MHz，25ms
    [bf_mat, angles] = fractionalBeamScan(pingData, fs, fc, interpFactor, angle_range);

    % 进行脉冲压缩
    bf_compressed = pulseCompress(bf_mat, fs*interpFactor, lfm_params);
    %% ---------- Upward‑facing polar‑style plot ----------
    newFS = 1000;
    dynRange = 20;  % dB
    env = abs(hilbert(bf_compressed));
    env = abs(downsampleBeamData(env, fs*interpFactor, newFS));
    BFdB = 20*log10(env./max(env(:))+eps);
    nSamples = size(BFdB,1);
    timeVec = single((0:nSamples-1)/newFS);  % 秒
    rangeVec = single(timeVec*c_water/2);              % 单程距离 (m)
    [AZG, RG] = meshgrid(single(angle_range), rangeVec);
    % 将方位角偏移 +90° 使 0° 指向 +Y 方向 → 扇面朝上
    [X,Y] = pol2cart(deg2rad(AZG + 90), RG);

    figure('Name','Fan Intensity Map – Upwards');
    surf(X,Y,BFdB,'EdgeColor','none');
    view(2); axis equal tight; shading interp;
    colormap(parula); caxis([-dynRange 0]);
    colorbar;
    xlabel('X (m)'); ylabel('Y (m)');

    title(sprintf('Upward DAS Fan Map – %d beams, %d samples',numel(angle_range),nSamples));
    pause(0.1)
end