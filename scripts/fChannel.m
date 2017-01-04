% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models the channel effects in the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% paths (Mx1 Integers) = Number of paths for each source in the system.
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths=[1;3;2]
% symbolsIn (MxR Complex) = Signals being transmitted in the channel
% delay (Cx1 Integers) = Delay for each path in the system starting with
% source 1
% beta (Cx1 Integers) = Fading Coefficient for each path in the system
% starting with source 1
% DOA = Direction of Arrival for each source in the system in the form
% [Azimuth, Elevation]
% SNR = Signal to Noise Ratio in dB
% array = Array locations in half unit wavelength. If no array then should
% be [0,0,0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% symbolsOut (FxN Complex) = F channel symbol chips received from each antenna
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ symbolsOut ] = fChannel( paths, symbolsIn, delay, beta, DOA, SNR_db, array )

numPaths = sum(paths);
numSources = size(symbolsIn, 1);
symbolsLength = size(symbolsIn, 2);
transmissionLength = symbolsLength + max(delay);
numSensors = size(array,1);

%% Calculate AWGN

rng('default'); % Reset rng
SNR_linear = 10^(SNR_db/10); % Find linear SNR
E_symbol = sum(abs(symbolsIn(1,:).^2))/symbolsLength; % Find symbol energy
N_0 = E_symbol/SNR_linear; % Find noise spectral density

noiseSigma = sqrt(N_0/2); % noise amplitude
awgn = noiseSigma*(randn(transmissionLength,numSensors) + 1j*randn(transmissionLength,numSensors));

%% Simulate Channel Paths

% zero fill where delays cause signal to end/start before/after others
symbolsAccumulated = zeros(transmissionLength,numSensors);

for pathIndex = 1:numPaths
    
    %determine which source this path comes from
    for srcIndex = 1:numSources
        firstPath = sum(paths(1:srcIndex-1)) + 1;
        lastPath = sum(paths(1:srcIndex));
        
        if pathIndex >= firstPath && pathIndex <= lastPath
            pathSrc = srcIndex;
            break
        end
    end
    
    symbolsPath_initial = symbolsIn(pathSrc,:).';    
    
    arrayManifoldVector = spv(array, DOA(pathIndex, :));
    
    % Add channel effects to symbols
    symbolsPath = (beta(pathIndex) .* symbolsPath_initial) * arrayManifoldVector.' ;
    
    % Add symbols to accumulated channel symbols
    symbolsStart = delay(pathIndex) + 1;
    symbolsEnd   = delay(pathIndex) + symbolsLength;
    symbolsAccumulated(symbolsStart:symbolsEnd,:) = symbolsAccumulated(symbolsStart:symbolsEnd,:) + symbolsPath;
    
end

%% Add noise to accumulated symbols

symbolsOut = symbolsAccumulated + awgn;

end