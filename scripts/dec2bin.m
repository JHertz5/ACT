function [ bin ] = dec2bin( dec, length)
% Efficient decimal integer to binary array  converter
% MSB first
% ~10000 times faster than de2bi
dec = double(dec); % convert to double, in case dec is an int
binLength = ceil(log2(dec+1));
bin = zeros(1,length,'uint8');

if binLength > length
    fprintf('ERROR: length %i is too small to contain binary vector of %i\n',length,dec)
end

for i = 1:binLength
    
    quotient = floor(dec/2);
    remainder = rem(dec,2);
    
    bin(length + 1 - i) = remainder;
    dec = quotient;
end

end

