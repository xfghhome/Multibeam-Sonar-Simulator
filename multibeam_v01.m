%% Init

clear
close all
clc
addpath(genpath('SourceCode'))

%% ---------- 适配水下多波束声纳参数 ----------
c_water = 1500;                 % 水下声速 (m/s)
f_c     = 900e3;                % 中心频率 900 kHz
BW      = 100e3;                % 调频带宽 100 kHz (850-950 kHz)
pulseT  = 5e-3;                 % 脉宽 5 ms
fs      = 2e6;                  % 2 MHz 采样率 (> 2×最高频)

%%   Load the surface and prepare the surface normals

nFreqsSim = 21;
structMeshPreparation = struct();
structMeshPreparation.orientation = [0 0 0];
structMeshPreparation.position = [0 0 0];
structMeshPreparation.vertexScaling = 1/1000;
structMeshPreparation.FLIPNORMALS = 1;
structMeshPreparation.vecFreqSim = linspace(f_c-BW/2, f_c+BW/2, nFreqsSim);
structMeshPreparation.fileNameMesh = "Data/Models/wall.stl";

structMeshPreparation.BRDFTransitionPosition = 0.4;
structMeshPreparation.BRDFTransitionSlope = 2;
structMeshPreparation.BRDFExponentSpecular = linspace(8,5,nFreqsSim);
structMeshPreparation.BRDFExponentDiffractive = linspace(70,70,nFreqsSim);

structMeshPreparation.materialStrengthSpecular = 10*linspace(1,0.8,nFreqsSim);
structMeshPreparation.materialStrengthDiffractive = 0.05*linspace(0.5,1,nFreqsSim);
structMeshPreparation.materialSTransitionPosition = structMeshPreparation.BRDFTransitionPosition;
structMeshPreparation.materialSTransitionSlope = 2;

structMeshPreparation.precomputeCurvature = 1;

structSurface = prepareMeshSurface(structMeshPreparation, 1);

%% Setup the structs for processing
%% ---------- 生成 128 阵元水平线阵 ----------

lambda   = c_water / f_c;      % 波长 ≈ 1.667 mm
d_elem   = lambda/2;           % λ/2 间距 ≈ 0.833 mm
N_elem   = 128;                % 阵元数
yCoords  = ((0:N_elem-1) - (N_elem-1)/2) * d_elem;   % 对称分布

% Make the receiver arrays
coordsReceivers = [ zeros(N_elem,1) , yCoords.' , zeros(N_elem,1) ]; % [x y z]
coordsEmitter = [ 0 0 0 ];

% Combine it all into the sensor
structSensor = struct();
structSensor.position        = [2 0 0].';   % 阵列参考点 (可根据场景修改)
structSensor.orientation     = [0 0 0];     % 水平朝 +X 方向
structSensor.coordsEmitter   = [ 1.98 0 0 ]; % 发射器位于阵前 2 cm 处
structSensor.coordsReceivers = coordsReceivers;
structSensor.nMics           = N_elem;

%% ---------- 更新仿真全局参数 ----------
structSimulationParameters                = struct();
structSimulationParameters.doPlot         = 0;  % 禁用所有绘图
structSimulationParameters.numSamplesImpresp = round(0.03*fs); % 30 ms 时窗
structSimulationParameters.sampleRateImpresp = fs;
structSimulationParameters.limitsAzimuth     = [-60 60];
structSimulationParameters.limitsElevation   = [-10 10];
structSimulationParameters.numberOfDirections        = 1e6;
structSimulationParameters.numberOfDirectionsPerCall = 2e5;
structSimulationParameters.vecFreqSim        = structMeshPreparation.vecFreqSim;
structSimulationParameters.numSamplesIRFilter= 256;
structSimulationParameters.IRFilterGaussAlpha= 5;
structSimulationParameters.numDiffractionPoints     = 0.5e5;
structSimulationParameters.approximateImpulseResponseCutDB = -90;
structSimulationParameters.approximateImpulseResponse      = 0;
structSimulationParameters.ditherRaytracing  = 1;
structSimulationParameters.speedOfSound      = c_water;

%% --------- 计算脉冲响应 ---------

tic
structSimulationResult = calculateImpulseResponseFast(structSensor, structSurface, structSimulationParameters);
toc

%% --------- 生成发射信号并合成通道数据 ---------

fStart = 900e3 - BW/2;          % 850 kHz
sigEmit = createLFM(BW, pulseT, fs, fStart, "up");
nSamples  = size(structSimulationResult.impulseResponse,1);

dataMics  = zeros(nSamples, N_elem);
for ch = 1:N_elem
    dataMics(:,ch) = conv(structSimulationResult.impulseResponse(:,ch), sigEmit, 'same');
end

%% ---------- Delay‑and‑Sum 波束形成 ----------

anglesDeg = -60:0.5:60;           % 波束扫描角
nAng      = numel(anglesDeg);
bfData    = zeros(nSamples, nAng);

% 预计算时间向量，便于插值实现分数延迟
sampleIdx = (0:nSamples-1).';
for a = 1:nAng
    theta = deg2rad(anglesDeg(a));
    % 各阵元延迟 (s)
    tau   = (yCoords.' * sin(theta)) / c_water;  % 1×128
    tauSamp = tau * fs;                          % 转为采样点
    for ch = 1:N_elem
        sig = dataMics(:,ch);
        bfData(:,a) = bfData(:,a) + interp1(sampleIdx, sig, sampleIdx - tauSamp(ch), 'linear', 0);
    end
end
bfData = bfData / N_elem;  % 归一化

% 包络检波以便显示
envBf  = abs(hilbert(bfData));

%% ---------- 绘图：角度‑距离强度图 ----------

rangeVec = sampleIdx * c_water / (2*fs); % 单程距离 (m)
figure;
imagesc(anglesDeg, rangeVec*1000, 20*log10(envBf / max(envBf(:))));
set(gca,'YDir','normal');
colormap(jet);
caxis([-60 0]);
colorbar;
xlabel('Azimuth (°)');
ylabel('Range (mm)');
title('Delay‑and‑Sum Beamformed Image');

grid on;
