clc
close all

%add data directory to path
if contains(pwd, 'ACT')
    dataPath = strcat( extractBefore(pwd, 'ACT'), 'ACT/data');
    addpath(char(dataPath));
else
    dataPath = ''; %dataPath is empty vector
    fprintf('Move to ACT directory\n');
end


%% Initialise Values

X = 8;  % H => 8
Y = 10; % J => 10

carrierFreq = 1;
carrierAmp = 1;
phi = X + 2*Y;

%% Generate PN sequences

pnPolynomialDegree = 4;
pn1Register = ones(pnPolynomialDegree, 1, 'logical');
pn2Register = ones(pnPolynomialDegree, 1, 'logical');

codeLength = 2^(pnPolynomialDegree) - 1;
pn1Seq = zeros(codeLength, 1);
pn2Seq = zeros(codeLength, 1);

for seqIndex = 1:codeLength

    pn1Seq(seqIndex) = pn1Register(4);
    pn2Seq(seqIndex) = pn2Register(4);
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

%% Generate Gold sequences

delayGold = 1 + mod(X + Y, 12); % all values above initial delayGold satisfy the inequality
goldSeq1 = fGoldSeq(pn1Seq, pn2Seq, delayGold);

while sum(goldSeq1, 1) ~= 8 % while gold code is not balanced
    fprintf('delayGold = %i does not provide balanced gold sequence\n', delayGold);
    delayGold = delayGold + 1;
    goldSeq1 = fGoldSeq(pn1Seq, pn2Seq, delayGold);
end

fprintf('delayGold = %i  is the smallest integer to satsfy inequality and provide balanced gold sequence\n ', delayGold);

goldSeq2 = fGoldSeq(pn1Seq, pn2Seq, delayGold + 1); % gold sequence for user 2 uses d + 1
goldSeq3 = fGoldSeq(pn1Seq, pn2Seq, delayGold + 2); % gold sequence for user 3 uses d + 2

%% Plot sequences

figure

subplot (2,2,1) % top left
stairs(pn1Seq, 'LineWidth', 2, 'Marker', 'o')
ylabel('PN Sequence 1')

subplot (2,2,3) % bottom left
stairs(pn2Seq, 'LineWidth', 2, 'Marker', 'o')
ylabel('PN Sequence 2')
xlabel('Sequence index')

subplot (3,2,2) % top right
stairs(goldSeq1, 'LineWidth', 2, 'Marker', 'o')
ylabel('Gold Seq 1')

subplot (3,2,4) % middle right
stairs(goldSeq2, 'LineWidth', 2, 'Marker', 'o')
ylabel('Gold Seq 2')

subplot (3,2,6) % bottom right
stairs(goldSeq3, 'LineWidth', 2, 'Marker', 'o')
ylabel('Gold Seq 3')
xlabel('Sequence index')

%% Read image into vector

image1 = imread('Photo1.jpg');
figure
imshow(image1);

%%

% if ~isempty(dataPath)
%     save(char(strcat(dataPath, '/wine_separatedData')),'training_classes','validation_classes','testing_classes','training_raw','validation_raw','testing_raw','training_norm','validation_norm','testing_norm')
% else
%     save('wine_separatedData','training_classes','validation_classes','testing_classes','training_raw','validation_raw','testing_raw','training_norm','validation_norm','testing_norm')
% end