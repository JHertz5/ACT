clc
close all
clearvars -EXCEPT showPlots SNR_dB

%add data directory to path
addpath('data')

fprintf('Q3\n')

%% Initialise Values

% SNR_dB
if ~exist('SNR_dB', 'var')
    SNR_dB = input('INPUT REQUIRED: Enter channel SNR level in dB: ');
    if isempty(SNR_dB)
        SNR_dB = 40;
        fprintf('\tNo value entered, SNR set to default value\n')
    end
end
fprintf('SNR set to %i dB\n', SNR_dB);

% showPlots
if ~exist('showPlots', 'var') || isempty(showPlots)
    showPlots = input('INPUT REQUIRED: Enter showPlots value: ');
    if isempty(SNR_dB)
        showPlots = 1;
        fprintf('\tshowPlots set to default value of 1\n')
    end
end
if showPlots == 0
    fprintf('showPlots == 0, no plots will be shown\n')
else
    fprintf('showPlots == %i\n', showPlots)
end

% phi
X = 8;  % H => 8
Y = 10; % J => 10
phi = X + 2*Y;
fprintf('X = %i, Y = %i -> phi = %i\n', X, Y, phi)

%% Read & Convert Photos

fprintf('Reading and converting images ... \n')

imgSize = 160*112*3*8;

figure

fprintf('\tReading image 1\n')
subplot(1,3,1)
[img1, imgW, imgH] = fImageSource('Photo1.jpg', imgSize);
fprintf('\tReading image 2\n')
subplot(1,3,2)
[img2, ~, ~] = fImageSource('Photo2.jpg', imgSize);
title('Original Images - 1, 2, 3')
fprintf('\tReading image 3\n')
subplot(1,3,3)
[img3, ~, ~] = fImageSource('Photo3.jpg', imgSize);
fprintf('\tComplete\n')

%% Transmitter Start

fprintf('Transmitter\n')

%% Generate M sequences

fprintf('\tGenerating M sequences ...\n')

MSeq1 = fMSeqGen([1 0 0 1 1]); % D^4 + D^1 + 1
MSeq2 = fMSeqGen([1 1 0 0 1]); % D^4 + D^3 + 1

fprintf('\t\tComplete\n')

%% Generate Gold sequences

fprintf('\tGenerating Gold sequences ...\n')

delayGold = 1 + mod(X + Y, 12); % all values above initial delayGold satisfy the inequality
goldSeq1 = fGoldSeq(MSeq1, MSeq2, delayGold);

while sum(goldSeq1, 1) ~= 8 % while gold code is not balanced
    fprintf('\t\tdelayGold = %i does not provide balanced gold sequence\n', delayGold);
    delayGold = delayGold + 1;
    goldSeq1 = fGoldSeq(MSeq1, MSeq2, delayGold);
end

fprintf('\t\tdelayGold = %i  is the smallest integer to satsfy inequality and provide balanced gold sequence\n', delayGold);

goldSeq2 = fGoldSeq(MSeq1, MSeq2, delayGold + 1); % gold sequence for user 2 uses d + 1
goldSeq3 = fGoldSeq(MSeq1, MSeq2, delayGold + 2); % gold sequence for user 3 uses d + 2

fprintf('\t\tComplete\n')

%% Perform DS-QPSK modulation 

fprintf('\tPerforming DS-QPSK modulation ...\n')

fprintf('\t\tModulating image 1\n')
symbols1 = fDSQPSKModulator(img1, goldSeq1, phi);
fprintf('\t\tModulating image 2\n')
symbols2 = fDSQPSKModulator(img1, goldSeq2, phi);
fprintf('\t\tModulating image 3\n')
symbols3 = fDSQPSKModulator(img1, goldSeq3, phi);

%% Plotting Transmitter Data

if showPlots >= 1
    
    figure
    scatter(real(symbols1), imag(symbols1), 75, 'x')
    line(xlim, [0 0], 'Color', 'black') % x axis line
    line([0 0], ylim, 'Color', 'black') % y axis line
    title('Constellation Diagram of QPSK symbols - Transmitter 1 Output')
    ylabel('Quadrature')
    xlabel('In Phase')
    grid on
end

fprintf('\t\tComplete\n')

%% Channel Start

fprintf('Channel\n')

%% Initialise Channel Values

fprintf('\tInitialising channel values ...\n')

paths = [ 1; 1; 1 ];
symbolsIn = (cat(2, symbols1, symbols2, symbols3))';
delay = [ 3; 7; 12 ];
beta  = [ 0.4; 0.7; 0.2 ];
DOA = [ 30 0; 90 0; 150 0 ];

% Set up sensor position vectors
arrayAngles = deg2rad(30:(360/5):318)'; % 5 sensors uniformly, circularly spaced, starting at 30 deg anticlockwise
% radius is selected such that the minimum inter-antenna spacing is 1
arrayRadius = sqrt(1/( (cos(arrayAngles(1))-cos(arrayAngles(2)))^2 + (sin(arrayAngles(1))-sin(arrayAngles(2)))^2 ));
array = [ arrayRadius*cos(arrayAngles) arrayRadius*sin(arrayAngles) zeros(5,1) ];

% Plot diagram of sensor positions
if showPlots >= 2
    
    figure
    hold on
    % Plot a radius circle
    t = linspace(0,2*pi,100);
    plot3(arrayRadius.*sin(t),arrayRadius.*cos(t),zeros(1,100),'Color','red');
    
    % Plot sensor points
    scatter3(array(:,1),array(:,2),array(:,3),'bx','LineWidth',2)
    text(array(:,1),array(:,2),array(:,3),{' Sensor 1',' Sensor 2',' Sensor 3',' Sensor 4',' Sensor 5'})

    line(xlim,[0 0],[0 0],'Color','black')
    line([0 0],ylim,[0 0],'Color','black')
    line([0 0],[0 0],zlim,'Color','black')
    title('Diagram of Sensor Positions')
    xlabel('x (\lambda/2)')
    ylabel('y (\lambda/2)')
    zlabel('z (\lambda/2)')
    view(-12,64)
end

fprintf('\t\tComplete\n')

%% Apply Channel Effects

fprintf('\tApplying channel effects ...\n')

symbolsChannelOut = fChannel(paths, symbolsIn, delay, beta, DOA, SNR_dB, array);

fprintf('\t\tComplete\n')

%% Receiver Start

fprintf('Receiver\n')

%% Receiver Detection
% Estimate the number of incoming signals using AIC criterion

fprintf('\tPerforming receiver detection ...\n')

symbolsOut = symbolsChannelOut.'; % get symbolsOut into correct orientation for processing

[numAntennas,numSamples] = size(symbolsOut);
covarianceMatrix = symbolsOut*symbolsOut'/numSamples; % find covariance mantrix
[eigenVecs,eigenVals] = eig(covarianceMatrix); %eigenvalue decomposition
eigenVals = sort(diag(eigenVals),'descend'); %sort the eigenvalue in the decreasing order

criterionAIC = zeros(1,numAntennas);
criterionMDL = zeros(1,numAntennas);

%Compute criterion
for antennaIndex = 0:numAntennas-1
    coef = 1/(numAntennas-antennaIndex);
    a = coef*sum(eigenVals(antennaIndex+1:numAntennas));
    g = prod(eigenVals(antennaIndex+1:numAntennas)).^(coef); 
    criterionAIC(antennaIndex+1) = -log(((g/a)^(numSamples*(numAntennas-antennaIndex))))+antennaIndex*(2*numAntennas-antennaIndex); 
    criterionMDL(antennaIndex+1) = -log(((g/a)^(numSamples*(numAntennas-antennaIndex))))+0.5*antennaIndex*(2*numAntennas-antennaIndex)*log(numSamples);
end    

[~,numSourcesEst] = min(criterionMDL);
numSourcesEst = numSourcesEst-1; % decrement 1 because the value antennaIndex=0 is stored in the index 1.

fprintf('\t\tNumber of incoming signals estimated as %i\n', numSourcesEst)
fprintf('\t\tComplete\n')

%% Receiver Estimation of DOA
% Estimate DOA of each signal using MUSIC algorithm

fprintf('\tPerforming channel parameter estimation ...\n')

eigVecsNoise = eigenVecs(:,1:end-numSourcesEst); % Noise eigenvectors

azRange = (0:1:180)'; % Azimuth values to search
elRange = (0:1:180)'; % Elevation values to search

% MUSIC spectrum computation
musicAmplitude = zeros(length(azRange),length(elRange));
for elIndex = 1:length(elRange)
    % Elevation search value
    el = elRange(elIndex);
    
    intersectionSearch = spv(array, [ azRange ones(size(azRange))*el ]);

    % Compute azimuth spectrum for this elevation
    musicAmplitude(:,elIndex) = sum(abs(intersectionSearch'*eigVecsNoise).^2,2);
end

musicAmplitude_dB = -10*log10(musicAmplitude.'/numAntennas);
 
if showPlots >= 1
    % Plot MUSIC spectrum
    figure();
    surf(azRange, elRange, musicAmplitude_dB);
    shading interp;
    title('MUSIC Spectrum');
    xlabel('Azimuth (degrees)');
    ylabel('Elevation (degrees)');
    zlabel('MUSIC spectrum (dB)');
    grid on; axis tight;
end

% find top 3 values to be DOA estimates
musicAmplitudeSorted = sort( reshape(musicAmplitude_dB,length(azRange)*length(elRange),1), 'descend');
musicTop3Threshold = musicAmplitudeSorted(4);

[ elMaxIndex, azMaxIndex ] = find(musicAmplitude_dB > musicTop3Threshold);
DOAest = [ azRange(azMaxIndex) elRange(elMaxIndex) ];
% DOAest(:,2) = [ 0; 0; 0 ];

fprintf('\t\tDOAs estimated as: \n')
for i = 1:size(DOAest,1)
    fprintf('\t\t\tSig %i: Az = %i deg, El = %i deg\n', i, DOAest(i,1), DOAest(i,2))
end

%% Receiver Estimation of Delay

testOffset = 0;
maxDelay = length(symbolsOut) - length(symbolsIn); % This is the greatest possible value of the delay

%  Would be finding Maximum Likelihood Estimator or MUSIC value for each
%  value of delay between 0 and maxDelay, and then taking the optimum
%  value, but I have not been able to figure out how to calculate the array
%  vector manifold for path delay

delayEst = delay(1);
fprintf('Note: Currently there is no delay estimation\n')

%% Receiver Reception
% Use weiner-hopf beamformer to receive signals

% Calculate Weiner-Hopf beamformer weights
% Would ideally have used superresolution beamformer, but not enough time
beamformWeightWH = inv(covarianceMatrix)*(spv(array,DOAest(1,:)));

% conjugate is needed to get correct answer, no idea why
symbolsReceived = conj(beamformWeightWH'*symbolsOut);

if showPlots >= 1
    % Plot beamformer information as a polar and euclidean plot
    az = 0:360;
    gain = zeros(size(az)); % only plotting in azimuth
    for azIndex = 1:length(az)
        gain(azIndex) = abs(beamformWeightWH'*spv(array,[az(azIndex) 0]));
    end
    gain = 10*log10(gain);
    
    gain_polarPlot = gain;
    gain_polarPlot(gain < 0) = 0;
    
    figure
    subplot(1,2,1)
    plot(az,gain)
    line([30 30],ylim)
    line([90 90],ylim)
    line([150 150],ylim)
    subplot(1,2,2)
    polar(deg2rad(az),gain_polarPlot)
end
%% Extracting Desired Signal

img1SymbolsStart = delayEst + 1;
img1SymbolsEnd   = length(symbols1) + delayEst;
img1SymbolsReceived = symbolsReceived(img1SymbolsStart:img1SymbolsEnd);

%% Plot Receiver input

if showPlots >= 4
    % Plot constellation diagram
    figure
    hold on
    scatter(real(symbolsReceived(:)), imag(symbolsReceived(:)), 'b.')
    line(xlim, [0 0], 'Color', 'black') % x axis line
    line([0 0], ylim, 'Color', 'black') % y axis line
    title('Constellation Diagram of QPSK symbols - Receiver Input at Sensor 1')
    ylabel('Imaginary Axis')
    xlabel('Real Axis')
    grid on
end

%% Perform DS-QPSK demodulation

fprintf('\tPerforming DS-QPSK de-modulation ...\n')
fprintf('\t\tDe-modulating image 1\n')
imgReceived1 = fDSQPSKDemodulator(img1SymbolsReceived, goldSeq1, phi);

fprintf('\t\tComplete\n')

%% Analysis of Receiver Data

errorRatio_img1 = length(find(img1 ~= imgReceived1))/length(img1);
fprintf('\t%f%% of received image values are incorrect (Channel SNR = %idB)\n',errorRatio_img1*100,SNR_dB)

if showPlots >= 1
    % Display images
    fprintf('\tPlotting received image 1\n')
    
    figure
    subplot(1,2,1)
    fImageSink(img1, imgSize, imgW, imgH)
    title('Original')
    subplot(1,2,2)
    fImageSink(imgReceived1, imgSize, imgW, imgH)
    title('Received')
    
end

%% End of Simulation

fprintf('Simulation complete\n')