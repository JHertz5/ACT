% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs channel estimation for the desired source using the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (Fx1 Complex) = R channel symbol chips received
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired source used in the modulation process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_estimate = Vector of estimates of the delays of each path of the
% desired signal
% DOA_estimate = Estimates of the azimuth and elevation of each path of the
% desired signal
% beta_estimate = Estimates of the fading coefficients of each path of the
% desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delay_estimate, DOA_estimate, beta_estimate]=fChannelEstimation(symbolsIn,goldSeq)

%% Discretiser and Manifold Extender

fprintf('\tPerforming spatiotemporal data transformation...\n')

seqLength = length(goldSeq);
[ numSamples, numAntennas ] = size(symbolsIn);

% Generate spaciotemporal samples

extensionLength = 2*seqLength;
X = zeros(numAntennas*extensionLength, numSamples);
tappedDelays = zeros(numAntennas, extensionLength);

% Fill tappedDelays
% for sampleIndex = 1:extensionLength
%     for tapIndex = 1:extensionLength-1
%         tappedDelays(:,tapIndex) = tappedDelays(:,tapIndex+1);
%     end
%     tappedDelays(:,end) = symbolsIn(sampleIndex,:);
% end

for sampleIndex = 1:numSamples    
    % step tapped delays
    for tapIndex = 1:extensionLength-1
        tappedDelays(:,tapIndex) = tappedDelays(:,tapIndex+1);
    end
    tappedDelays(:,end) = symbolsIn(sampleIndex,:);
    X(:,sampleIndex) = reshape(tappedDelays.',[],1); % store vectorised data
end

fprintf('\t\tComplete\n')

%% Spatiotemporal Channel Estimator

fprintf('\tPerforming channel parameter estimation ...\n')

covarianceMatrix = X*X'/numSamples; % find covariance mantrix
[eigenVecs,eigenVals] = eig(covarianceMatrix); %eigenvalue decomposition
eigenVals = diag(eigenVals);
plot(eigenVals)

noiseSubspaceRange = 1 : find(eigenVals<0.02,1,'last');
eigVecsNoise = eigenVecs(:,noiseSubspaceRange); % Noise eigenvectors

toaRange = (0:1:seqLength-1)'; % Elevation values to search

packedPNcode = [goldSeq;zeros(size(goldSeq))];
% Generate shifting matrix
shiftMatrix = zeros(2*seqLength);
shiftMatrix(2:end,1:end-1) = eye(2*seqLength-1);

% MUSIC spectrum computation
musicAmplitude = zeros(size(toaRange));
for toaIndex = 1:length(toaRange)
    % Elevation search value
    toa = toaRange(toaIndex);
    
%     spaceManifoldVector = spv(array, [ doaRange zeros(size(doaRange)) ]);
%     spaceManifoldVector = 1;
    timeManifoldVector  = (shiftMatrix^toa)*packedPNcode;
%     intersectionSearch  = kronProd(spaceManifoldVector,timeManifoldVector);
%     intersectionSearch  = kronProd(spaceManifoldVector,timeManifoldVector); % Kronecker Product
%     intersectionSearch = spaceManifoldVector * timeManifoldVector; % Kronecker Product
    
    %kronProd
    %kron

    % Compute azimuth spectrum for this elevation
    musicAmplitude(toaIndex) = sum(abs(timeManifoldVector'*eigVecsNoise).^2,2);
end

musicAmplitude_dB = -10*log10(musicAmplitude.'/numAntennas);
 
% if showPlots >= 1
    % Plot MUSIC spectrum
    figure;
    plot(toaRange, musicAmplitude_dB);
    shading interp;
    title('MUSIC Spectrum');
    xlabel('Azimuth (degrees)');
%     ylabel('TOA');
    ylabel('MUSIC spectrum (dB)');
    grid on; axis tight;
% end

% find top 3 values to be DOA estimates
% musicAmplitudeSorted = sort( reshape(musicAmplitude_dB,length(doaRange)*length(toaRange),1), 'descend');
% musicTop3Threshold = musicAmplitudeSorted(4);

% [ toaMaxIndex, doaMaxIndex ] = find(musicAmplitude_dB > musicTop3Threshold);
% DOA_estimate = doaRange(doaMaxIndex) 
% delay_estimate = toaRange(toaMaxIndex);

% fprintf('\t\tDOAs estimated as: \n')
% for i = 1:size(DOA_estimate,1)
%     fprintf('\t\t\tSig %i: Az = %i deg, El = %i deg\n', i, DOA_estimate(i), delay_estimate(i))
% end

end