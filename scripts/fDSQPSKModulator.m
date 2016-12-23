% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform DS-QPSK Modulation on a vector of bits using a gold sequence
% with channel symbols set by a phase phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% bitsIn (Px1 Integers) = P bits of 1's and 0's to be modulated
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence to be used in the modulation process
% phi (Integer) = Angle index in degrees of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% symbolsOut (Rx1 Complex) = R channel symbol chips after DS-QPSK Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ symbolsOut ] = fDSQPSKModulator( bitsIn, goldSeq, phi )

% chipPeriod and carrierPeriod are given nominal values
T_carrier = 5; % carrier period
T_chip = 1; % chip period
chipsPerBit = T_carrier/T_chip; % chipsPerBit must be an integer value

symbolsLength = length(bitsIn)/2;
chipsLength = symbolsLength * chipsPerBit;
seqLength = length(goldSeq);

symbolsPaired = zeros(symbolsLength, 1, 'double');
symbolsModulated = complex( symbolsPaired, symbolsPaired ); % initialise as a symbolsBinary sized complex vector

%% Perform QPSK modulation

phi_rad = deg2rad(phi);

for symbolIndex = 1:symbolsLength
    
    symbolsPaired(symbolIndex) = bitsIn(symbolIndex*2 - 1) * 2 + bitsIn(symbolIndex*2); % group bits into pairs
    
    %check for invalid values
    if symbolsPaired(symbolIndex) < 0 || symbolsPaired(symbolIndex) > 3
        fprintf('ERROR: bitsInPaired contained value outside range 0-3: %i @ index %i\n', bitsInPaired(symbolIndex), symbolIndex)
        return
    end
    
    % Find symbols in polar form
    theta = phi_rad + symbolsPaired(symbolIndex)*pi/2; % NOTE: diagram was vague on what phi should be for n ~= 00
    r = sqrt(2);
    
    % Convert to a + bi form
    symbolsModulated(symbolIndex) = complex(r*cos(theta), r*sin(theta));
end

%% Spread using Gold Sequence
% As specified in SomeNotes.1, "The messages are first modulated using QPSK,
% then spread by the gold sequences"

% Change gold sequence form
goldSeq = 1 - 2*goldSeq; % 1 -> -1, 0 -> 1

symbolsOut = complex(zeros(chipsLength, 1), zeros(chipsLength, 1));

for chipIndex = 1:chipsLength
    seqIndex = mod(chipIndex-1, seqLength) + 1; % index for gold sequence
    symbolIndex = ceil(chipIndex/chipsPerBit); % index for QPSK symbols
    
    symbolsOut(chipIndex) = goldSeq(seqIndex) * symbolsModulated(symbolIndex);
end

end