
%% Init

    clear
    close all
    clc
    addpath(genpath('SourceCode'))
%% ---------- 适配水下多波束声纳参数 ----------
    c_water = 1500;                 % 水下声速 (m/s)
    f_c     = 900e3;               % 中心频率 900 kHz
    BW      = 100e3;               % 调频带宽 100 kHz (850-950 kHz)
    pulseT  = 5e-3;                % 脉宽 5 ms
    fs      = 2e6;                 % 4 MHz 采样率 ( > 2×最高频 )
%%   Load the surface and prepare the surface normals
    
    nFreqsSim = 21;
    structMeshPreparation = struct();
    structMeshPreparation.orientation = [ 0 0 0];
    structMeshPreparation.position = [ 0 0 0];
    structMeshPreparation.vertexScaling = 1/1000;
    structMeshPreparation.FLIPNORMALS = 1;
    structMeshPreparation.vecFreqSim = linspace(f_c-BW/2, f_c+BW/2, nFreqsSim); 
    structMeshPreparation.fileNameMesh = "Data/Models/wall.stl";

    structMeshPreparation.BRDFTransitionPosition = 0.4;
    structMeshPreparation.BRDFTransitionSlope = 2;
    structMeshPreparation.BRDFExponentSpecular = linspace( 8,5,nFreqsSim);
    structMeshPreparation.BRDFExponentDiffractive = linspace( 70,70,nFreqsSim);

    structMeshPreparation.materialStrengthSpecular = 10*linspace( 1,0.8,nFreqsSim);
    structMeshPreparation.materialStrengthDiffractive = 0.05*linspace( 0.5,1,nFreqsSim);
    structMeshPreparation.materialSTransitionPosition = structMeshPreparation.BRDFTransitionPosition;
    structMeshPreparation.materialSTransitionSlope = 2;
    
    structMeshPreparation.precomputeCurvature = 1;
    
    structSurface = prepareMeshSurface( structMeshPreparation, 1 );

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
    structSimulationParameters.doPlot         = 0;
    structSimulationParameters.numSamplesImpresp = round(0.03*fs); % 30 ms 时窗
    structSimulationParameters.sampleRateImpresp = fs;
    structSimulationParameters.limitsAzimuth     = [-60 60];  % 声纳波束扇区更宽
    structSimulationParameters.limitsElevation   = [-10 10];  % 基本水平面
    structSimulationParameters.numberOfDirections        = 10e5;
    structSimulationParameters.numberOfDirectionsPerCall = 2e5;
    structSimulationParameters.vecFreqSim        = structMeshPreparation.vecFreqSim;
    structSimulationParameters.numSamplesIRFilter= 256;
    structSimulationParameters.IRFilterGaussAlpha= 5;
    structSimulationParameters.numDiffractionPoints     = 0.5e5;
    structSimulationParameters.approximateImpulseResponseCutDB = -90;
    structSimulationParameters.approximateImpulseResponse      = 0;
    structSimulationParameters.ditherRaytracing  = 1;
    structSimulationParameters.speedOfSound      = c_water;

%% Plot the setup, to make sure normals are OK

    nNormals = size( structSurface.surfaceNormals, 1 );
    normalSpacing = round( nNormals / 1000);
    idxNormalsPlot = 1 : normalSpacing : nNormals;
    figure(1212)       
        clf
        hp = patch('faces', structSurface.surfaceFaces, 'vertices', structSurface.surfaceVertices, 'FaceColor', [ 1 0 0], 'EdgeAlpha', 0.4); 
        xlabel( 'X-Axis' )
        ylabel( 'Y-Axis' )
        zlabel( 'Z-Axis' )
        grid on
        
        hold on;
            drawTriad( structSensor.position(:), structSensor.orientation(:), 0.05)
        hold off;

        hold on
        for cntNormal = 1 : length( idxNormalsPlot )
            curIdxPlot = idxNormalsPlot( cntNormal );
            curPointOrigin = structSurface.surfaceVertices( curIdxPlot, : );
            curPointNormal = curPointOrigin + structSurface.surfaceNormals( curIdxPlot, : ) * 0.02;
            curPointsPlot = [ curPointOrigin ; curPointNormal ];
            plot3( curPointsPlot(:,1), curPointsPlot(:,2), curPointsPlot(:,3), 'g' );
        end
        hold off;
        axis equal
        set( gca, 'view', [21.6026   19.6988]);

%% Process the raytracing
  
    tic
    structSimulationResult = calculateImpulseResponseFast( structSensor, structSurface, structSimulationParameters );
    toc

%% Process the impulseresponses  

    sigEmit = fm_sweep( 150e3, 40e3, 450e3, 0.4, 1, 10 );
    dataMics = zeros( size( structSimulationResult.impulseResponse ) );
    for cntMic = 1 : size(dataMics,2)
        dataMics( :, cntMic ) = conv( structSimulationResult.impulseResponse(:, cntMic), sigEmit, 'same' );
    end

%% Plot the final result, including the reflection points.

    strengthsReflections = structSimulationResult.strengthsReflected; 
    pointsReflections = structSimulationResult.pointsReflected;
    nPointsReflected = size( pointsReflections, 1 );
    plotPointsSkip = 10;

    idxFrequencyPoints = 1;
    idxMicPoints = 1;
    strengthsReflectionsNorm = normLin( squeeze( strengthsReflections( :, idxFrequencyPoints, idxMicPoints ) ) ); 

    dbCutoffPlot = 60;
    strengthsReflectionsNorm = normLog( strengthsReflectionsNorm, - dbCutoffPlot );
    strengthsReflectionsNorm = normLin( strengthsReflectionsNorm + dbCutoffPlot );
    
    strengthsReflectionsIndexer = round( strengthsReflectionsNorm *254 ) +1;
    cmapPlot = parula(256);
    
    colorsPoints = cmapPlot( strengthsReflectionsIndexer, : );
        figure()       
        
            subplot( 2,3,[1 2 4 5])
                cla
                hp = patch('faces', structSurface.surfaceFaces, 'vertices', structSurface.surfaceVertices, 'FaceColor', [ 0.2 0.3 0.3 ], 'EdgeAlpha', 0.3); 
                xlabel( 'X-Axis' )
                ylabel( 'Y-Axis' )
                zlabel( 'Z-Axis' )
                grid on
                set( gca, 'view', [67.1313 7.3250] );
                axis equal

                hold on;
                    drawTriad( structSensor.position(:), structSensor.orientation(:), 0.01)
                    for cntPoint = 1 : plotPointsSkip: nPointsReflected
                        plot3( pointsReflections(cntPoint, 1), pointsReflections(cntPoint,2), pointsReflections(cntPoint,3), '.', 'Color', colorsPoints(cntPoint,:), 'MarkerSize', strengthsReflectionsNorm(cntPoint)*20+0.1 )
                    end
                hold off;

                % Create a colorbar
                c = colorbar;
                c.Label.String = 'Normalized Reflection Strength based on BRDF (dB)';
                colormap(cmapPlot);
                caxis([-dbCutoffPlot, 0]);
                colorbarPosition = [0.1, 0.1, 0.02, 0.8];
                set(c, 'Position', colorbarPosition);
                c.Label.FontSize = 14;  % Adjust the font size as needed
                c.FontSize = 12;
                set( gcf, 'position', [273         348        1725         990])
                title('Reflectivity plot of scene with sensor position/orientation', 'fontsize', 16)

                camlight
                hp.FaceLighting = 'gouraud';
                hp.AmbientStrength = 0.4;
                hp.DiffuseStrength = 0.8;
                hp.SpecularStrength = 0.1;
                hp.SpecularExponent = 1;
                hp.BackFaceLighting = 'unlit';
                
            
            numSamplesPlotIR = 15000;
            IRLeft =  structSimulationResult.impulseResponse(1:numSamplesPlotIR,1);
            IRRight =  structSimulationResult.impulseResponse(1:numSamplesPlotIR,2);
            timeVec = (1:numSamplesPlotIR)/structSimulationParameters.sampleRateImpresp;
            subplot(3,3,3)
                plot( timeVec*1000,IRLeft, 'linewidth', 1.5 )
                grid on
                xlabel( 'Time (ms)')
                ylabel( 'Impulse Response (pressure - au)')
                title( 'Left Impulse Response' )
            subplot(3,3,6)
                plot( timeVec*1000,IRRight, 'linewidth', 1.5 )                
                grid on
                xlabel( 'Time (ms)')
                ylabel( 'Impulse Response (pressure - au)')
                title( 'Right Impulse Response' )
                
            subplot(3,3,9)
                plot( timeVec*1000,IRRight, 'linewidth', 1.5 ) 
                hold on
                plot( timeVec*1000,IRLeft, 'linewidth', 1.5 )
                hold off;
                grid on
                xlabel( 'Time (ms)')
                ylabel( 'Impulse Response (pressure - au)')
                title( 'Impulse Response Overlay' )              
                legend( 'Right', 'Left')

