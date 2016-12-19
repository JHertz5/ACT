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

%% Read Photo

fprintf('Extracting\n')
img1_raw = imread('Photo1.jpg');
img2_raw = imread('Photo2.jpg');
img3_raw = imread('Photo3.jpg');

figure
subplot(1,3,1)
imshow(img1_raw)
subplot(1,3,2)
imshow(img2_raw)
subplot(1,3,3)
imshow(img3_raw)

imgHeight = size(img1_raw, 1);
imgWidth  = size(img1_raw, 2);
imgDepth  = size(img1_raw, 3);
layerSize = imgWidth * imgHeight;

img1_vector = zeros(imgHeight*imgWidth*imgDepth, 1, 'uint8');
img1_reconstructed = zeros(imgHeight, imgWidth, imgDepth, 'uint8');
img2_vector = zeros(imgHeight*imgWidth*imgDepth, 1, 'uint8');
img2_reconstructed = zeros(imgHeight, imgWidth, imgDepth, 'uint8');
img3_vector = zeros(imgHeight*imgWidth*imgDepth, 1, 'uint8');
img3_reconstructed = zeros(imgHeight, imgWidth, imgDepth, 'uint8');

fprintf('Vectorising\n')
for k = 1:imgDepth
    for j = 1:imgWidth
        colStart = (k-1)*layerSize + (imgHeight)*(j-1)+1;
        colEnd   = (k-1)*layerSize + (imgHeight)*(j);
        
        img1_vector(colStart:colEnd) = img1_raw(:,j,k);
        img2_vector(colStart:colEnd) = img2_raw(:,j,k);
        img3_vector(colStart:colEnd) = img3_raw(:,j,k);
    end
end

vectorLength = size(img1_vector, 1);

img1_binVector = zeros(8*vectorLength, 1);
img1_binReconstructed = zeros(8*vectorLength, 1);
img2_binVector = zeros(8*vectorLength, 1);
img2_binReconstructed = zeros(8*vectorLength, 1);
img3_binVector = zeros(8*vectorLength, 1);
img3_binReconstructed = zeros(8*vectorLength, 1);

fprintf('Converting decimal to binary\n')
for i = 1:vectorLength
    byteStart = (i-1)*8 + 1;
    byteEnd   = i*8;
    img1_binVector(byteStart:byteEnd) = de2bi( img1_vector(i),8 );
    img2_binVector(byteStart:byteEnd) = de2bi( img2_vector(i),8 );
    img3_binVector(byteStart:byteEnd) = de2bi( img3_vector(i),8 );
end

fprintf('Converting binary to decimal\n')
for i = 1:vectorLength
    byteStart = (i-1)*8 + 1;
    byteEnd   = i*8;
    img1_binReconstructed(i) = bi2de( img1_binVector(byteStart:byteEnd)' );
    img2_binReconstructed(i) = bi2de( img2_binVector(byteStart:byteEnd)' );
    img3_binReconstructed(i) = bi2de( img3_binVector(byteStart:byteEnd)' );
end

fprintf('Reconstructing matrix\n')
for k = 1:imgDepth
    for j = 1 : imgWidth
        colStart = (k-1)*layerSize + (imgHeight)*(j-1)+1;
        colEnd   = (k-1)*layerSize + (imgHeight)*(j);
        img1_reconstructed(:,j,k) = img1_binReconstructed(colStart:colEnd);
        img2_reconstructed(:,j,k) = img2_binReconstructed(colStart:colEnd);
        img3_reconstructed(:,j,k) = img3_binReconstructed(colStart:colEnd);
    end
end

if img1_raw ~= img1_reconstructed
    fprintf('Image 1 reconstruction is incorrect');
end

if img2_raw ~= img2_reconstructed
    fprintf('Image 2 reconstruction is incorrect');
end

if img3_raw ~= img3_reconstructed
    fprintf('Image 3 reconstruction is incorrect');
end

figure
subplot(1,3,1)
imshow(img1_reconstructed)
subplot(1,3,2)
imshow(img2_reconstructed)
subplot(1,3,3)
imshow(img3_reconstructed)

