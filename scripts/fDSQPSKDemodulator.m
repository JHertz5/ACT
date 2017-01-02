% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform demodulation of the received data using <INSERT TYPE OF RECEIVER>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (Fx1 Integers) = R channel symbol chips received
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired signal to be used in the demodulation process
% phi (Integer) = Angle index in degrees of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bitsOut (Px1 Integers) = P demodulated bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ bitsOut ] = fDSQPSKDemodulator( symbolsIn, goldSeq, phi )

seqLength = length(goldSeq); % also the number of chips per bit
symbolsLength = length(symbolsIn);
bitPairsLength = symbolsLength/seqLength;

bitPairsModulated = complex(zeros(bitPairsLength, 1), zeros(bitPairsLength, 1));
bitsOut = zeros(bitPairsLength*2, 1, 'uint8');

%% De-Spread using Gold Sequence

% Change gold sequence form
goldSeq = 1 - 2*goldSeq; % 1 -> -1, 0 -> 1

for symbolIndex = 1:symbolsLength
    seqIndex = mod(symbolIndex-1, seqLength) + 1; % index for gold sequence
    bitPairsIndex = ceil(symbolIndex/seqLength); % index for QPSK symbols
    
    % Accumulate chips to be averaged
    bitPairsModulated(bitPairsIndex) = bitPairsModulated(bitPairsIndex) + goldSeq(seqIndex) * symbolsIn(symbolIndex);    
    
    % Every full gold sequence cycle, take the mean of the accumulation
    if mod(symbolIndex, seqLength) == 0
        bitPairsModulated(bitPairsIndex) = bitPairsModulated(bitPairsIndex) / seqLength;
    end
end

%% Perform QPSK de-modulation

phi_rad = deg2rad(phi);

for bitPairsIndex = 1:bitPairsLength
    
    % Find phase angle
    theta = angle(bitPairsModulated(bitPairsIndex));
    if theta < 0
        theta = theta + 2*pi;
    end
    
    % Convert phase angle to binary values    
    % This is the decision device. Signals have equal weight and 
    % probability, so decision device is very simple
    bitPair = mod(round((theta - phi_rad)*2/pi),4);

    % Split bit pairs into bits
    bitsOut(2*bitPairsIndex - 1 : 2*bitPairsIndex) = dec2bin(bitPair, 2);
end

end
