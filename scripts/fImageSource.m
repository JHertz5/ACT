% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads an image file with AxB pixels and produces a column vector of bits
% of length Q=AxBx3x8 where 3 represents the R, G and B matrices used to
% represent the image and 8 represents an 8 bit integer. If P>Q then
% the vector is padded at the bottom with zeros.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% filename (String) = The file name of the image
% P (Integer) = Number of bits to produce at the output - Should be greater
% than or equal to Q=AxBx3x8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bitsOut (Px1 Integers) = P bits (1's and 0's) representing the image
% x (Integer) = Number of pixels in image in x dimension
% y (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bitsOut,imgW,imgH]=fImageSource(filename,P)

fprintf('fImageSource - \nfilename = %s\n', filename)

%% Extracting Image

img_raw = imread(filename);
imshow(img_raw)

imgH = size(img_raw, 1); % height
imgW = size(img_raw, 2); % width
imgD = size(img_raw, 3); % depth should be 3 for R G B

layerSize = imgW * imgH;
vectorLength = imgH*imgW*3*8; %Q

% Error checking
if imgD ~= 3
    fprintf('ERROR: Image depth is %i, but should be 3\n', imgD)
    return
elseif P < vectorLength % P < Q
    fprintf('ERROR: Image size is too large for specified number of output bits %i \n%i is the correct number of bits\n', P, vectorLength);
    return
end

%% Vectorising Image

bitsOut = zeros(P, 1, 'uint8'); % using P rather than Q to implement zero padding

% column-wise raster scan of the image
for k = 1 : imgD
    for j = 1 : imgW
        for i = 1 : imgH
            byteStart = ((k-1)*layerSize + (imgH)*(j-1) + (i-1))*8 + 1;
            byteEnd   = ((k-1)*layerSize + (imgH)*(j-1) + i)*8;

            bitsOut(byteStart:byteEnd) = de2bi( img_raw(i,j,k),8 );
        end
    end
end

fprintf('%s successfully read\n\n', filename);

end