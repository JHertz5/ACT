clc
close all

%1:17:50 timestap for panopto presentation in final lecture
% Direct Sequence QPSK

%add data directory to path
if contains(pwd, 'ACT')
    dataPath = strcat( extractBefore(pwd, 'ACT'), 'ACT/data');
    addpath(char(dataPath));
else
    dataPath = ''; %dataPath is empty vector
    fprintf('Move to ACT directory\n');
end

%% Create PN codes and find gold codes

pnPolynomialDegree = 4;
pn1Register = ones(pnPolynomialDegree, 1, 'logical');
pn2Register = ones(pnPolynomialDegree, 1, 'logical');

codeLength = 2^(pnPolynomialDegree) - 1;
pn1Code = zeros(codeLength, 1);
pn2Code = zeros(codeLength, 1);

for codeIndex = 1:codeLength

    pn1Code(codeIndex) = pn1Register(4);
    pn2Code(codeIndex) = pn2Register(4);
    pn1Input =  mod( (pn1Register(3) + pn1Register(4)), 2 ); % D^4 + D + 1
    pn2Input =  mod( (pn2Register(1) + pn2Register(4)), 2 ); % D^4 + D^3 + 1

    % working backwards, shift each data into next register
    for regIndex = pnPolynomialDegree:-1:2
        pn1Register(regIndex) = pn1Register(regIndex-1); 
        pn2Register(regIndex) = pn2Register(regIndex-1); 
    end

    pn1Register(1) = pn1Input;
    pn2Register(1) = pn2Input;
end

% if ~isempty(dataPath)
%     save(char(strcat(dataPath, '/wine_separatedData')),'training_classes','validation_classes','testing_classes','training_raw','validation_raw','testing_raw','training_norm','validation_norm','testing_norm')
% else
%     save('wine_separatedData','training_classes','validation_classes','testing_classes','training_raw','validation_raw','testing_raw','training_norm','validation_norm','testing_norm')
% end