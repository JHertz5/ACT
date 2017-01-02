function [ dec ] = bin2dec( bin )
% Efficient binary array to decimal integer converter
% ~64 times faster than bi2de

dec = 0;
for i = 1 : length(bin)
    dec = dec + bin(i) * 2^(length(bin) - i);
end
end

