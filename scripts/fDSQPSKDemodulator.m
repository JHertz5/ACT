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

% chipPeriod and carrierPeriod are given nominal values
T_carrier = 5; % carrier period
T_chip = 1; % chip period
chipsPerBit = T_carrier/T_chip; % chipsPerBit must be an integer value

chipsLength = length(symbolsIn);
symbolLength = chipsLength/chipsPerBit;
seqLength = length(goldSeq);

%% De-Spread using Gold Sequence

% Change gold sequence form
goldSeq = 1 - 2*goldSeq; % 1 -> -1, 0 -> 1

symbolsModulated = complex(zeros(symbolLength, 1), zeros(symbolLength, 1));

for chipIndex = 1:chipsLength
    seqIndex = mod(chipIndex-1, seqLength) + 1; % index for gold sequence
    symbolIndex = ceil(chipIndex/chipsPerBit); % index for QPSK symbols
    
    % Accumulate chips to be averaged later
    symbolsModulated(symbolIndex) = symbolsModulated(symbolIndex) + goldSeq(seqIndex) * symbolsIn(chipIndex);    
    
    if mod(chipIndex, chipsPerBit) == 0
        symbolsModulated(symbolIndex) = symbolsModulated(symbolIndex) / chipsPerBit;
    end
end

%% Perform QPSK de-modulation

phi_rad = deg2rad(phi);
bitsOut = zeros(symbolLength*2, 1, 'uint8');

for symbolIndex = 1:symbolLength
    
    % Find phase angle
    theta = angle(symbolsModulated(symbolIndex));
    if theta < 0
        theta = theta + 2*pi;
    end
    
    % Convert phase angle to binary values    
    bitPair = uint8((theta - phi_rad)*2/pi);


    % check for invalid values
    if bitPair < 0 || bitPair > 3
        fprintf('ERROR: bitsInPaired contained value outside range 0-3: %i @ index %i\n', symbolsPaired(symbolIndex), symbolIndex)
        return
    end
    
    % Split bit pairs into bits
    bitsOut(2*symbolIndex - 1 : 2*symbolIndex) = de2bi(bitPair, 2, 'left-msb');
end

end
