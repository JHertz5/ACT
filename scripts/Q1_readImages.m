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

%% Read & Convert Photos

imgSize = 430080;

figure
subplot(1,3,1)
[img1, imgW, imgH] = fImageSource('Photo1.jpg', imgSize);
subplot(1,3,2)
[img2, ~, ~] = fImageSource('Photo2.jpg', imgSize);
subplot(1,3,3)
[img3, ~, ~] = fImageSource('Photo3.jpg', imgSize);

figure
subplot(1,3,1)
fImageSink(img1, imgSize, imgW, imgH)
subplot(1,3,2)
fImageSink(img2, imgSize, imgW, imgH)
subplot(1,3,3)
fImageSink(img3, imgSize, imgW, imgH)

%% Save Data

if ~isempty(dataPath)
    save(char(strcat(dataPath, '/ImageVectors')),'img1','img2','img3','imgW','imgH','imgSize')
else
    save('ImageVectors','img1','img2','img3','imgW','imgH','imgSize')
end