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

bitPairsLength = length(bitsIn)/2;
seqLength = length(goldSeq); % also the number of chips per bit
symbolsLength = bitPairsLength * seqLength;

bitPairs = zeros(bitPairsLength, 1, 'double');
bitPairsModulated = complex( bitPairs, bitPairs ); % initialise as a symbolsBinary sized complex vector

%% Perform QPSK modulation

phi_rad = deg2rad(phi);

for bitPairsIndex = 1:bitPairsLength
    
    bitPairs(bitPairsIndex) = bitsIn(bitPairsIndex*2 - 1) * 2 + bitsIn(bitPairsIndex*2); % group bits into pairs with bin to dec conversion
    
    %check for invalid values
    if bitPairs(bitPairsIndex) < 0 || bitPairs(bitPairsIndex) > 3
        fprintf('ERROR: bitsInPaired contained value outside range 0-3: %i @ index %i\n', bitsInPaired(bitPairsIndex), bitPairsIndex)
        disp('TODO: replace with errror function')
        return
    end
    
    % Find symbols in polar form
    theta = phi_rad + bitPairs(bitPairsIndex)*pi/2;
    r = sqrt(2);
    
    % Convert to a + bi form
    bitPairsModulated(bitPairsIndex) = complex(r*cos(theta), r*sin(theta));
end

%% Spread using Gold Sequence
% As specified in SomeNotes.1, "The messages are first modulated using QPSK,
% then spread by the gold sequences"

% Change gold sequence form
goldSeq = 1 - 2*goldSeq; % 1 -> -1, 0 -> 1

symbolsOut = complex(zeros(symbolsLength, 1), zeros(symbolsLength, 1));

for symbolIndex = 1:symbolsLength
    seqIndex = mod(symbolIndex-1, seqLength) + 1; % index for gold sequence
    bitPairsIndex = ceil(symbolIndex/seqLength); % index for QPSK symbols
    
    symbolsOut(symbolIndex) = goldSeq(seqIndex) * bitPairsModulated(bitPairsIndex);
end

end