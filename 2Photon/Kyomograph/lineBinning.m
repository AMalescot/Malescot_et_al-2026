function [bKymo] = lineBinning(kymo,bin)

nLine = length(kymo);
x = size(kymo,2);
nBin = floor(nLine/bin);

bKymo = zeros(nBin,x,'uint16');
for i = 1:nBin
    bKymo(i,:) = mean(kymo(1+(i-1)*bin:i*bin,:),1);
end

