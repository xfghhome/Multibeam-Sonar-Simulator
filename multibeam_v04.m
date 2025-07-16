%% Multi‑beam sonar simulation with DAS beamforming – fan display (upwards)
% -------------------------------------------------------------
% (c) 2025 Kyume – generates impulse responses via ray‑tracing,
% applies 128‑channel Delay‑and‑Sum beamforming, draws an upward‑
% facing fan (polar‑style) intensity map, and visualises the scene
% geometry with the sonar array elements marked as red dots.
% -------------------------------------------------------------

clear; close all; clc;
addpath(genpath('SourceCode'));

%% ---------- System parameters ----------
c_water = 1500;              % Speed of sound in water (m/s)
f_c     = 900e3;             % Centre frequency (Hz)
eps = 0.1;
BW      = 100e3;             % Sweep bandwidth (Hz)
pulseT  = 5e-3;              % Pulse length (s)
fs      = 2e6;               % Sampling rate (Hz)
interpFactor = 2;
angle_range = -60:1:60;    % 扫描角度范围
% simulation frequencies
nFreqsSim = 3;

%% ---------- Geometry / STL surface ----------
structMeshPreparation = struct();
structMeshPreparation.orientation            = [0 0 180];
structMeshPreparation.position               = [10 3 0];
structMeshPreparation.vertexScaling          = 1/100;  % mm→m
structMeshPreparation.FLIPNORMALS            = 1;
structMeshPreparation.vecFreqSim             = linspace(f_c-BW/2,f_c+BW/2,nFreqsSim);
structMeshPreparation.fileNameMesh           = "Data/Models/ball.stl";
structMeshPreparation.BRDFTransitionPosition = 0.4;
structMeshPreparation.BRDFTransitionSlope    = 2;
structMeshPreparation.BRDFExponentSpecular   = linspace(8,5,nFreqsSim);
structMeshPreparation.BRDFExponentDiffractive= linspace(70,70,nFreqsSim);
structMeshPreparation.materialStrengthSpecular   = 10*linspace(1,0.8,nFreqsSim);
structMeshPreparation.materialStrengthDiffractive= 0.05*linspace(0.5,1,nFreqsSim);
structMeshPreparation.materialSTransitionPosition= structMeshPreparation.BRDFTransitionPosition;
structMeshPreparation.materialSTransitionSlope   = 2;
structMeshPreparation.precomputeCurvature    = 1;

structSurface = prepareMeshSurface(structMeshPreparation,0);

%% ---------- Linear 128‑element array ----------
lambda  = c_water/f_c;
d_elem  = lambda/2;                    % half‑wavelength spacing
N_elem  = 128;
yCoords = ((0:N_elem-1)-(N_elem-1)/2)*d_elem; % centred
coordsReceivers = [zeros(N_elem,1) , yCoords.' , zeros(N_elem,1)];

structSensor = struct();
structSensor.position        = [0 0 0].';    % array reference (m)
structSensor.orientation     = [0 0 0];
structSensor.coordsEmitter   = [0.02 0 0];    % emitter 2 cm in front (absolute)
structSensor.coordsReceivers = coordsReceivers; % relative to position
structSensor.nMics           = N_elem;

%% ---------- Simulation parameters ----------
structSimulationParameters = struct();
structSimulationParameters.doPlot           = 0;
structSimulationParameters.numSamplesImpresp= round(0.025*fs); % 25 ms
structSimulationParameters.sampleRateImpresp= fs;
structSimulationParameters.limitsAzimuth    = [-80 80];
structSimulationParameters.limitsElevation  = [-80 80];
structSimulationParameters.numberOfDirections        = 1e6;
structSimulationParameters.numberOfDirectionsPerCall = 2e5;
structSimulationParameters.vecFreqSim       = structMeshPreparation.vecFreqSim;
structSimulationParameters.numSamplesIRFilter       = 256;
structSimulationParameters.IRFilterGaussAlpha       = 5;
structSimulationParameters.numDiffractionPoints     = 5e4;
structSimulationParameters.approximateImpulseResponseCutDB = -90;
structSimulationParameters.approximateImpulseResponse      = 0;
structSimulationParameters.ditherRaytracing  = 1;
structSimulationParameters.speedOfSound      = c_water;

%% ---------- Calculate impulse responses ----------
structSimulationResult = calculateImpulseResponseFast(structSensor,structSurface,structSimulationParameters);

%% ---------- Generate LFM pulse & convolve ----------
fStart = f_c - BW/2;                               % 850 kHz
sigEmit = createLFM(BW,pulseT,fs,fStart,"up");   % 实值 LFM

dataMics = zeros(size(structSimulationResult.impulseResponse),"single");
for m = 1:N_elem
    dataMics(:,m) = conv(structSimulationResult.impulseResponse(:,m),sigEmit,'same');
end


%% ---------- Scene visualisation (geometry + array) ----------
figure('Name','Scene Geometry');
% Surface mesh
hp = patch('Faces',structSurface.surfaceFaces,'Vertices',structSurface.surfaceVertices, ...
           'FaceColor',[0.2 0.3 0.3],'EdgeAlpha',0.3);
hold on;
% Sonar array elements (absolute positions)
arrayPos = structSensor.position.' + coordsReceivers; % Nx3
scatter3(arrayPos(:,1),arrayPos(:,2),arrayPos(:,3),20,'r','filled');
% Emitter (yellow)
scatter3(structSensor.coordsEmitter(1),structSensor.coordsEmitter(2),structSensor.coordsEmitter(3),50,'y','filled');
% Coordinate triad for sensor orientation (if drawTriad exists)
% if exist('drawTriad','file')
%     drawTriad(structSensor.position(:),structSensor.orientation(:),0.05);
% end
hold off;
axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
view([45 20]);
camlight; hp.FaceLighting='gouraud'; hp.AmbientStrength=0.4; hp.DiffuseStrength=0.8;

title('Simulation Scene: surface & sonar array');

%% ---------- Upward‑facing polar‑style plot ----------
[bf_mat, angles] = fractionalBeamScan(dataMics, fs, f_c, interpFactor, angle_range);

% 进行脉冲压缩
lfm_params.f0 = f_c-BW/2;        % 起始频率
lfm_params.f1 = f_c+BW/2;        % 终止频率
lfm_params.T_p = pulseT;        % 脉冲宽度 5 ms
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


