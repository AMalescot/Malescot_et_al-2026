function [av_S1]  = AverageTrials(S,fps,final_tBin, chanTag, TrialType)
%% Description: 
% function to average single trial of the same type, and create a .tiff
% video
%% Input:
% - S: pre-processed signal
% - fps: framerate
% - final_tBin: final temporal binning
% - chanTag: tag of the channel
% - TrialType: 'Real'
%% Output:
% - av_S1: post-binning averaged trials
%
% Written by Éric Martineau - Université de Montreal
%% Mean of all trials by type %%
av_S1 = mean(S(:,:,:,TrialType),4,"omitnan");

if sum(TrialType)<size(S,4)
    av_S1(:,:,:,2) = mean(S(:,:,:,~TrialType),4,"omitnan");
end
clear S

%% Temporal averaging %%
binSize = fps / final_tBin;
if binSize < 1 % if final_tBin > fps, don't bin 
    binSize = 1;
else
    binSize = round(binSize); %if final_tBin < fps, bin frame to get as close as possible to desired fps
end

A = zeros(size(av_S1,1),size(av_S1,2),floor(size(av_S1,3)/binSize),size(av_S1,4));

for i = 1:size(A,3)
    A(:,:,i,:) = mean(av_S1(:,:,(binSize*(i-1)+1):i*binSize,:),3,"omitnan");      
end
clear C D

%% Save images for Fiji %%
binTag = int2str(fps/binSize);

saveFloatToTiff_ioi(A(:,:,:,1),strcat(chanTag, "_Real_", binTag, "hzBin"));
if size(A,4)>1
    saveFloatToTiff_ioi(A(:,:,:,2),strcat(chanTag, "_Sham_", binTag, "hzBin"));
end