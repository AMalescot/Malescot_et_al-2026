function [Signal,Mask_donut] = Analysis_donuts_timbitsApp(AllSignalsFilt,Parameters,options)
%% Descrition: 
%function to extract the signal from a disk and several surrounding donuts. 
% You have to go through Step1 and 2 before to use that function.
%
%% Input: - AllSignals : struct get with IOI Analysis
%        - Parameters: struct get with IOI Analysis with all the settings
%        - options : struct with number of iterations, radius of
%        disk/donut and threshold to segment the pick of neuronal activity
%% Output: - Signal: struct (Real and Sham) with number of iterations as row
%         - Mask_donut: xy-nIterations matrix of the donut masks.
%
% Written by Antoine Malescot, adapted by Éric Martineau - Université de Montreal
%%

nb_of_donut = options.iterations;
radius = options.radius; %radius in pixels, ok if not integer because calculates the pixels number over the value.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Donuts surrrounding the central activity %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dims = size(Parameters.windowMask);
if size(Parameters.Mask,3) == 1 %if only one signal based ROI, take this one
    [row, col] = ind2sub(dims,find(Parameters.Mask == 1)); %find all pixels in signal-based ROI
    Centroid = [mean(col) mean(row)]; %xy center of activity

    [columnsInImage, rowsInImage] = meshgrid(1:dims(2), 1:dims(1));
    tmp_mask = (rowsInImage- Centroid(2)).^2 ...
        + (columnsInImage - Centroid(1)).^2 <= radius.^2; % Draw a circular mask centered on the centroid

    Mask_donut(:,:,1) = tmp_mask;
end

fields = fieldnames(AllSignalsFilt);
iterations = size(AllSignalsFilt,1)*nb_of_donut; %full length of waitbar

for j = 1:nb_of_donut
    if j>1
        %%% Create a circular mask with several surrounding donuts %%%
        [columnsInImage, rowsInImage] = meshgrid(1:dims(2), 1:dims(1));
        tmp_mask = (rowsInImage- Centroid(2)).^2 ...
            + (columnsInImage - Centroid(1)).^2 <= ((j)*radius).^2; % Draw a circular mask % Removed the +15 from antoine on the initial diameter
        Mask_donut(:,:,j) = logical(tmp_mask-max(Mask_donut(:,:,1:j-1),[],3)); % Create mask, or the donuts %Antoine had sum() here, although should be fine, max() is safer.
    elseif j == 1
        %find the center of neuronal analysis
        mask = Parameters.Mask;
        [row, col] = ind2sub(dims,find(mask == 1)); %find all pixels in signal-based ROI
        Centroid = [mean(col) mean(row)]; %xy center of activity

        [columnsInImage, rowsInImage] = meshgrid(1:dims(2), 1:dims(1));
        tmp_mask = (rowsInImage- Centroid(2)).^2 ...
            + (columnsInImage - Centroid(1)).^2 <= radius.^2; % Draw a circular mask centered on the centroid

        Mask_donut(:,:,1) = tmp_mask;
    end

    %%% Extract and assign the signals to a structure %%%
    for k = 1:length(fields)
        if ~isempty(AllSignalsFilt(1).(fields{k}))
            % Real
            tmp = sum(Mask_donut(:,:,j).*AllSignalsFilt(1).(fields{k}).*Parameters.windowMask,[1 2])./sum(Mask_donut(:,:,j),'all');
            Signal.Real.(fields{k})(j,:) = movmean(reshape(tmp,1,[]),[Parameters.RunAv-1 0]);
        end
        if ~isempty(AllSignalsFilt(2).(fields{k}))
            % Sham
            tmp = sum(Mask_donut(:,:,j).*AllSignalsFilt(2).(fields{k}).*Parameters.windowMask,[1 2])./sum(Mask_donut(:,:,j),'all');
            Signal.Sham.(fields{k})(j,:) = movmean(reshape(tmp,1,[]),[Parameters.RunAv-1 0]);
        end
    end
end
Mask_donut = Mask_donut.*Parameters.windowMask;
end

