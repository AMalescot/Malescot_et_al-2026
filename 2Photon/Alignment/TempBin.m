function [bstack] = TempBin(stack, bin)
% function binning frame across time to increase SNR
%
%% Input: 
%  - stack : time series, must be in XYTC format
%  - bin : number of frames to bin
%
%% Output:
%  - bstack: temporal binned data
%
% Written by Éric Martineau - Université de Montréal

%% Extract size
dims = size(stack);

bstack = zeros(dims(1), dims(2), floor(dims(3)/bin),size(stack,4));

%% Loop
for c = 1:size(stack,4)  
    for i = 1:floor(dims(3)/bin)
        idx = (i-1)*bin+1;
        bstack(:,:,i,c) = mean(stack(:,:,idx:idx+(bin-1),c),3);
    end
end

bstack = uint16(bstack);
