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

function [ symbolsOut ] = fChannel( paths, symbolsIn, delay, beta, DOA, SNR, array )

numPaths = sum(paths);
numSources = size(symbolsIn, 1);
symbolsLength = length(symbolsIn);

transmissionLength = symbolsLength + max(delay);
% zero fill where delays cause signal to end/start before/after others
symbolsAccumulated = zeros(1,transmissionLength);

for srcIndex = 1:numSources
    
    symbolsChannel = zeros(numPaths, symbolsLength);
    
    for pathIndex = 1:paths(srcIndex) 
        symbolsChannel(pathIndex,:) = symbolsIn(srcIndex,:);
    end
    
    % Find paths of source
    firstPath = sum(paths(1:srcIndex-1)) + 1;
    lastPath = sum(paths(1:srcIndex));
    srcPaths = (firstPath : lastPath)'; % used to get the delay, beta, and DOA of paths

    arrayManifoldVector = spv(array, DOA(srcPaths, :));
    S_H = (conj(arrayManifoldVector))';
    
    for pathIndex = 1:paths(srcIndex)
        srcPathIndex = srcPaths(pathIndex);
        symbolsChannel(pathIndex,:) = symbolsChannel(pathIndex,:) * S_H * beta(srcPathIndex);
        
        symbolsStart = delay(srcPathIndex);
        symbolsEnd   = delay(srcPathIndex) + symbolsLength;
        symbolsAccumulated(symbolsStart:symbolsEnd) = symbolsChannel(pathIndex,:);
    end
    
%     symbolsOut = symbolsAccumulated + 
    
end




%% Delay


%% DOA Adjustment


%% Fading


%% MAI



end