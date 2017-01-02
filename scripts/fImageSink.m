% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the received image by converting bits back into R, B and G
% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% bitsIn (Px1 Integers) = P demodulated bits of 1's and 0's
% Q (Integer) = Number of bits in the image
% x (Integer) = Number of pixels in image in x dimension
% y (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fImageSink(bitsIn,~,imgW,imgH)

imgD = 3; % R G B = 3 layers
img_reconstruction = zeros(imgH, imgW, imgD, 'uint8');

layerSize = imgW * imgH;

for k = 1 : imgD
    for j = 1 : imgW
        for i = 1 : imgH
           byteStart = ((k-1)*layerSize + (imgH)*(j-1) + (i-1))*8 + 1;
            byteEnd   = ((k-1)*layerSize + (imgH)*(j-1) + i)*8;
            img_reconstruction(i,j,k) = bin2dec( bitsIn(byteStart:byteEnd)' );
        end
    end
end

imshow(img_reconstruction)

end