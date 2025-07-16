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
BW      = 100e3;             % Sweep bandwidth (Hz)
pulseT  = 5e-3;              % Pulse length (s)
fs      = 2e6;               % Sampling rate (Hz)
% simulation frequencies
nFreqsSim = 3;
meshx = 10;
meshy = linspace(-3,3,11);
meshz = 0;
for i = 1:length(meshy)

    %% ---------- Geometry / STL surface ----------
    structMeshPreparation = struct();
    structMeshPreparation.orientation            = [0 0 180];
    structMeshPreparation.position               = [meshx meshy(i) meshz];
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
    savePingData(dataMics,['./OutputChannelData/ping_',num2str(i),'.dat'])

end