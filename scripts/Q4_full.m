clc
close all
clearvars -EXCEPT showPlots

%add data directory to path
addpath('data')

fprintf('Q4\nReceiver\n')

%% Initialise Values

% showPlots
if ~exist('showPlots', 'var') || isempty(showPlots)
    showPlots = input('INPUT REQUIRED: Enter showPlots value: ');
    if isempty(showPlots)
        showPlots = 1;
        fprintf('\tshowPlots set to default value of 1\n')
    end
end
if showPlots == 0
    fprintf('showPlots == 0, no plots will be shown\n')
else
    fprintf('showPlots == %i\n', showPlots)
end

% load personalised data
load jjh113.mat

% establish known quantities and parameters
numSources = length(Beta_1);
[ numAntennas, numSamples ] = size(Xmatrix);

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

%% Generate M sequences

fprintf('\tGenerating M sequences ...\n')

MSeq1 = fMSeqGen([1 0 0 1 0 1]); % D^5 + D^2 + 1
MSeq2 = fMSeqGen([1 1 1 1 0 1]); % D^5 + D^4 + D^3 + D^2 + 1

fprintf('\t\tComplete\n')

%% Generate Gold sequences

fprintf('\tGenerating Gold sequence ...\n')

goldSeq = fGoldSeq(MSeq1, MSeq2, phase_shift);

fprintf('\t\tComplete\n')

%% Discretiser and Manifold Extender

fprintf('\tPerforming spatiotemporal data transformation...\n')

seqLength = length(goldSeq);

% Generate spaciotemporal samples

extensionLength = 2*seqLength;
X = zeros(numAntennas*extensionLength, numSamples);
tappedDelays = zeros(numAntennas, extensionLength);

% Fill tappedDelays
for sampleIndex = 1:62
    for tapIndex = 1:extensionLength-1
        tappedDelays(:,tapIndex) = tappedDelays(:,tapIndex+1);
    end
    tappedDelays(:,end) = Xmatrix(:,sampleIndex);
end

for sampleIndex = 1:numSamples
    X(:,sampleIndex) = reshape(tappedDelays.',[],1); % store vectorised data
    
    % step tapped delays
    for tapIndex = 1:extensionLength-1
        tappedDelays(:,tapIndex) = tappedDelays(:,tapIndex+1);
    end
    tappedDelays(:,end) = Xmatrix(:,sampleIndex);
end

fprintf('\t\tComplete\n')

%% Spatiotemporal Channel Estimator

fprintf('\tPerforming channel parameter estimation ...\n')

covarianceMatrix = X*X'/numSamples; % find covariance mantrix
[eigenVecs,eigenVals] = eig(covarianceMatrix); %eigenvalue decomposition
eigenVals = diag(eigenVals);
% plot(eigenVals)

noiseSubspaceRange = 1 : find(eigenVals < 0.001,1,'last');
eigVecsNoise = eigenVecs(:,1:noiseSubspaceRange); % Noise eigenvectors

doaRange = (0:1:360)'; % Azimuth values to search
toaRange = (0:1:seqLength-1)'; % Elevation values to search

packedPNcode = [goldSeq;zeros(size(goldSeq))];
% Generate shifting matrix
shiftMatrix = zeros(2*seqLength);
shiftMatrix(2:end,1:end-1) = eye(2*seqLength-1);

% MUSIC spectrum computation
musicAmplitude = zeros(length(doaRange),length(toaRange));
for toaIndex = 1:length(toaRange)
    % TOA search value
    toaTest = toaRange(toaIndex);
    
    for doaIndex = 1:length(doaRange)
        doaTest = doaRange(doaIndex);
        
        spaceManifoldVector = spv(array, [ doaTest 0 ]);
        timeManifoldVector  = (shiftMatrix^toaTest)*packedPNcode;
        intersectionSearch  = kronProd(spaceManifoldVector,timeManifoldVector); % Kronecker Product
        
        % Compute azimuth spectrum for this elevation
        musicAmplitude(doaIndex,toaIndex) = sum(abs(intersectionSearch'*eigVecsNoise).^2,2);
    end
end

musicAmplitude_dB = -10*log10(musicAmplitude.'/numAntennas);
 
if showPlots >= 1
    % Plot MUSIC spectrum
    figure;
    surf(doaRange, toaRange, musicAmplitude_dB);
    shading interp;
    title('MUSIC Spectrum');
    xlabel('Azimuth (degrees)');
    ylabel('TOA');
    zlabel('MUSIC spectrum (dB)');
    grid on; axis tight;
end

% find top value to be DOA-TOA estimate
musicAmplitudeSorted = sort( reshape(musicAmplitude_dB,length(doaRange)*length(toaRange),1), 'descend');
musicTopThreshold = musicAmplitudeSorted(2);

[ toaMaxIndex, doaMaxIndex ] = find(musicAmplitude_dB > musicTopThreshold);
DOAest = doaRange(doaMaxIndex);
TOAest = toaRange(toaMaxIndex);

fprintf('\t\tEstimates: \n')
fprintf('\t\t\tDOA = %i deg, TOA = %i Tc\n', DOAest, TOAest)
fprintf('\t\tComplete\n')

%% Spatiotemporal Beamformer

fprintf('\tReceiving using beamformer ...\n')

spaceManifoldVector = spv(array, [ DOAest 0 ]);
timeManifoldVector  = (shiftMatrix^TOAest)*packedPNcode;
beamformWieghts  = kron(spaceManifoldVector,timeManifoldVector)*Beta_1(1); % Kronecker Product

Xreceived = (X.'*beamformWieghts)./Beta_1(1);
fprintf('\t\tComplete\n')

%% Demodulate
fprintf('\tPerforming DS-QPSK de-modulation ...\n')

msgLength = 8*60/2*seqLength;
msgRange = TOAest+1 : TOAest+msgLength;
msgReceived = fDSQPSKDemodulator(Xreceived(msgRange), goldSeq, phi_mod);
fprintf('\t\tComplete\n')

%% Decoding message
fprintf('\tDecoding message ...\n')

msgChars = [];
for charIndex = 1:60
    byteRange = 8*(charIndex-1)+1:8*charIndex;
    msgChars = [msgChars char(bin2dec(msgReceived(byteRange)))];
end

fprintf('\t\tMessage decoded as:\n\t\t%s\n', msgChars)
fprintf('\t\tComplete\n')

%% End of Simulation

fprintf('Simulation complete\n')