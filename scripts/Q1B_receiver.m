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

if ~(exist('showPlots', 'var') && showPlots == true)
    fprintf('showPlots is not true\n')
end

fprintf('Q1A - DS-QPSK Receiver\n')

%% Initialise Values

load ImageVectors
load Q1A_transmitterOutput

X = 8;  % H => 8
Y = 10; % J => 10
phi = X + 2*Y;
fprintf('\tX = %i, Y = %i -> phi = %i\n', X, Y, phi)

%% Generate M sequences

fprintf('\tGenerating M sequences ...\n')

MSeq1 = fMSeqGen([1 0 0 1 1]); % D^4 + D^1 + 1
MSeq2 = fMSeqGen([1 1 0 0 1]); % D^4 + D^3 + 1

fprintf('\t\tComplete\n')

%% Generate Gold sequences

fprintf('\tGenerating Gold sequences ...\n')

delayGold = 1 + mod(X + Y, 12); % all values above initial delayGold satisfy the inequality
goldSeq1 = fGoldSeq(MSeq1, MSeq2, delayGold);

while sum(goldSeq1, 1) ~= 8 % while gold code is not balanced
    fprintf('\t\tdelayGold = %i does not provide balanced gold sequence\n', delayGold);
    delayGold = delayGold + 1;
    goldSeq1 = fGoldSeq(MSeq1, MSeq2, delayGold);
end

fprintf('\t\tdelayGold = %i  is the smallest integer to satsfy inequality and provide balanced gold sequence\n', delayGold);

goldSeq2 = fGoldSeq(MSeq1, MSeq2, delayGold + 1); % gold sequence for user 2 uses d + 1
goldSeq3 = fGoldSeq(MSeq1, MSeq2, delayGold + 2); % gold sequence for user 3 uses d + 2

if exist('showPlots', 'var') && showPlots == true  
    figure
    
    subplot (2,2,1) % top left
    stairs(MSeq1, 'LineWidth', 2, 'Marker', 'x')
    ylabel('PN Sequence 1')
    xlim([1 15])
    
    subplot (2,2,3) % bottom left
    stairs(MSeq2, 'LineWidth', 2, 'Marker', 'x')
    ylabel('PN Sequence 2')
    xlabel('Sequence index')
    xlim([1 15])
    
    subplot (3,2,2) % top right
    stairs(goldSeq1, 'LineWidth', 2, 'Marker', 'x')
    ylabel('Gold Seq 1')
    xlim([1 15])
    
    subplot (3,2,4) % middle right
    stairs(goldSeq2, 'LineWidth', 2, 'Marker', 'x')
    ylabel('Gold Seq 2')
    xlim([1 15])
    
    subplot (3,2,6) % bottom right
    stairs(goldSeq3, 'LineWidth', 2, 'Marker', 'x')
    ylabel('Gold Seq 3')
    xlabel('Sequence index')
    xlim([1 15])
end

fprintf('\t\tComplete\n')

%% Perform DS-QPSK demodulation

fprintf('\tPerforming DS-QPSK de-modulation ...\n')
fprintf('\t\tDe-modulating image 1\n')
imgReceived1 = fDSQPSKDemodulator(symbols1, goldSeq1, phi);
fprintf('\t\tDe-modulating image 2\n')
imgReceived2 = fDSQPSKDemodulator(symbols2, goldSeq2, phi);
fprintf('\t\tDe-modulating image 3\n')
imgReceived3 = fDSQPSKDemodulator(symbols3, goldSeq3, phi);

if exist('showPlots', 'var') && showPlots == true
    % Plot constellation diagram
    figure
    plot( real(symbols1), imag(symbols1), 'x' )
    line([-1.5  1.5], [0 0], 'Color', 'black') % x axis line
    line([0 0], [-1.5  1.5], 'Color', 'black') % y axis line
    title('Constellation Diagram of QPSK symbols')
    ylabel('Imaginary Axis')
    xlabel('Real Axis')
    grid on
    
    % Display images
    figure
    subplot(1,3,1)
    fprintf('\t\tPlotting image 1\n')
    fImageSink(img1, imgSize, imgW, imgH)
    subplot(1,3,2)
    fprintf('\t\tPlotting image 2\n')
    fImageSink(img2, imgSize, imgW, imgH)
    subplot(1,3,3)
    fprintf('\t\tPlotting image 3\n')
    fImageSink(img3, imgSize, imgW, imgH)
end

fprintf('\t\tComplete\n')
