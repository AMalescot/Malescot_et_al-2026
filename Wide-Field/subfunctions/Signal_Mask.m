function [mask] = Signal_Mask(data, Threshold, Start, End)
%% Description:
% This function display a thresholded average signal.
%% Input:
% - data: data for mask extraction
% - Threshold: percentage of the normalized signal
% - Start: First frame for averaging
% - End: Last frame for averaging
%
%% Output:
% - mask: binary mask
%
% Written by Éric Martineau and Antoine Malescot - Université de Montréal

data = data(:,:,(Start+1):End);
data = mean(data,3,'omitnan');

maxVal = max(data,[],'all');
mask = data > (maxVal*Threshold);

%Prompt user to select right ROI
fig = figure();
if str2double(extractBetween(version,'(R',lettersPattern + ')')) > 2024
    theme('light');
end
imshow(mat2gray(~mask,[-1 1]));
[column, row] = ginput(1);
mask = bwselect(mask,column,row); 

clear max 
close(fig)


