% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes polynomial weights and produces an M-Sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% coeffs (Px1 Integers) = Polynomial coefficients. For example, if the
% polynomial is D^5+D^3+D^1+1 then the coeffs vector will be [1;0;1;0;1;1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% MSeq (Wx1 Integers) = W bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ MSeq ] = fMSeqGen( coeffs )

MPolynomialDegree = length(coeffs) - 1;
MRegister = ones(MPolynomialDegree, 1, 'uint8');

seqLength = 2^(MPolynomialDegree) - 1;
MSeq = zeros(seqLength, 1);

feedbackIndices = coeffs(2:end) == 1;
for seqIndex = 1:seqLength
    
    MSeq(seqIndex) = MRegister(MPolynomialDegree);
    
    MInput = mod(sum(MRegister(feedbackIndices)),2);
    
    % working backwards, shift each data into next register
    for regIndex = MPolynomialDegree:-1:2
        MRegister(regIndex) = MRegister(regIndex-1);
    end
    
    MRegister(1) = MInput;
end

end