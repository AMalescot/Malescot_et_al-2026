function [traces] = pullROIsignal(sig, Masks, avBin)
%% Description:
% function to extract the signal from inside the activity Mask
%
%% Input:
% - sig: structure containing preprocessed data
% - Masks: mask extracted from activity in Step2
% - avBin: number of frames to use for the moving average 
% - fps: framerate
%
%% Output:
% - traces: structure containing the extracted traces activity inside the
%   ROI signal
%
% Written by Éric Martineau and Antoine Malescot - Université de Montréal
%%
fields = fieldnames(sig);
Mask = logical(Masks);

for i = 1:length(fields)
    time = size(sig.(fields{i}),3);
    for t = 1:time
        frame = sig.(fields{i})(:,:,t);
        traces.(fields{i})(1,t) = mean(frame(Mask),'omitnan'); %series-roi-t matrix
    end
    traces.(fields{i}) = movmean(traces.(fields{i}), [(avBin-1) 0],'omitnan'); % Moving average 
end


    
