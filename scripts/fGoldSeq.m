% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes two M-Sequences of the same length and produces a gold sequence by
% adding a delay and performing modulo 2 addition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% mseq1 (Wx1 Integer) = First M-Sequence
% mseq2 (Wx1 Integer) = Second M-Sequence
% shift (Integer) = Number of chips to shift second M-Sequence to the right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% GoldSeq (Wx1 Integer) = W bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ goldSeq ] = fGoldSeq( MSeq1, MSeq2, shift )

seqLength = length(MSeq1); % W

goldSeq = zeros(seqLength, 1); % pre-allocate memory

for seqIndex = 1:seqLength
    % find delayed index 
    if seqIndex + shift > seqLength
        seq2Index = mod(seqIndex + shift, seqLength); % if statement to ensure that index 0 is skipped
    else
        seq2Index = seqIndex + shift;
    end
    
    goldSeq(seqIndex) = mod(MSeq1(seqIndex) + MSeq2(seq2Index), 2);
end

end