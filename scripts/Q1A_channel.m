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

fprintf('Q1A - DS-QPSK Channel\n')

%% Initialise Values

% Load pre-processed input images
load Q1A_transmitterOutput.mat

X = 8;  % H => 8
Y = 10; % J => 10
phi = X + 2*Y;
fprintf('\tX = %i, Y = %i -> phi = %i\n', X, Y, phi)

paths = [ 1; 1; 1 ];
symbolsIn = cat(2, symbols1, symbols2, symbols3);
delay = [ 3; 7; 12 ];
beta  = [ 0.4; 0.7; 0.2 ];
DOA = [ 30 0; 90 0; 150 0 ];
SNR = 0; %or 40
array = [ 0 0 0 ];

symbolsOut = fChannel(paths, symbolsIn, delay, beta, DOA, SNR, array);

%%

