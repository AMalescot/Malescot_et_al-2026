%% Note: 
% Function to extract the vascular diameter using Vasometrics. To use it,
% you need to be in your experiment folder and load your 2P time series
% data.

%% Input:
% AllData : structure containing all the video and metadata already
% preprocessed
% VasChan : channel number of vessel labeling (Alexa680)
% VesselID : text string to identify vessel
% Run batch process on one vessel segment at the time using different
% vesselID
%
% Written by Eric Martineau, adapted by Antoine Malescot - Université de Montreal

%% Inputs

VasChan = 2; % Channel of the vascular signal
VesselID = 'PA'; % Vessel ID
stimOnset = 4; % Start of the stimulation (in seconds)
penVessel = 0; % Measure a penetrating arteriole dilation (boolean value)

% Trial type label
RealID = "Real"; % Real
ShamID = "Sham"; % Sham

%% Create Vasometrics output folder
root = pwd;
mkdir Vasometrics;

out = strcat(root,'\Vasometrics\');

%% Batch vasometrics cmd and average table %%
close all
CLspace = 0.5; % Space between each cross lines (in um)
minDist = 0; % 0 = take the maximal FWHM, 1= take minimal FWHM
cLength = 100; % For standard method: input crossline length here if automatic not working. Empty array ([]) if you want auto.
cLengthR = 20; % For rosas method: length (in pixels) to add on each side of the diameter to determine crossline length.
alpha = 10; % Rotation angle between crossline for penetrating vessels, min value of 10o. 

% Preallocation
diam = zeros(size(AllData,1), size(AllData(1).Frames,3));

%Loop
for k = 1:length(AllData)
    % Run through Vasometrics
    disp(strcat("Analyzing trial",num2str(k))); 
    if k == 1
        if penVessel == 0 
            [diam(k,:), fwhms, CLs, fig] = Vasometrics(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, {}, CLspace, minDist, cLength,[]);
        elseif penVessel == 1 
            [diam(k,:), fwhms, CLs, fig] = Vasometrics_penetrating(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, {}, minDist, cLengthR,alpha,[]);
        end
        savefig(fig,fullfile(out,strcat(VesselID,'Crosslines')));
    else
        if penVessel == 0
            [diam(k,:), fwhms, ~, ~] = Vasometrics(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, CLs, CLspace, minDist, cLength,[]);
        elseif penVessel == 1
            [diam(k,:), fwhms, ~, ~] = Vasometrics_penetrating(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, CLs, minDist, cLengthR,alpha,[]);
        end
    end
    FWHMS{k} = fwhms;
end

%Store results
Parameters = struct('FinalFramerate', AllData(1).metadata.FrameRate,'CLspace',CLspace,'VesselChan',VasChan);
save(fullfile(out,strcat('diameters',VesselID,'.mat')),'diam','CLs','FWHMS','Parameters','-mat');

%% Plot diameter and save the figure

trialFlag = [AllData.Type]; % Trial type
fBin = round(0.05/ (1/AllData(1).metadata.FrameRate)); %number of frames in bin
nBin = floor(length(diam)/fBin); % number of bin
bDiam = zeros(size(diam,1),nBin); % preallocation of diameter value
t = zeros(1,nBin); % preallocation of time value

for i = 1:nBin % Bining of diameters 
    bDiam(:,i) = mean(diam(:,((i-1)*fBin+1):i*fBin),2); 
    t(1,i) = (i*fBin)/AllData(1).metadata.FrameRate; %timepoint = end of bin
end


fig = figure('Position', [256,128,1280,768]);
avDiam = zeros(length(bDiam),2);
relAvDiam = zeros(length(bDiam),2);

figure(fig)
Y1 = bDiam(trialFlag == "Real",:); % Real
Y2 = bDiam(trialFlag ~= "Real",:); % Sham
subplot(2,2,1);
plot(permute(t,[2,1]), permute(Y1,[2,1]));
title(strcat('Indiv Trials ',RealID));
xlabel('Time (s)');
ylabel('Diameter (um)');        

subplot(2,2,2);
if isempty(Y2) == 0
    plot(permute(t,[2,1]), permute(Y2,[2,1]));   
end
title(strcat('Indiv Trials ',ShamID));
xlabel('Time (s)');
ylabel('Diameter (um)');

subplot(2,2,3);
avDiam(:,1) = mean(Y1,1,'omitnan');
plot(permute(t,[2,1]), avDiam(:,1));
if isempty(Y2) == 0
    hold on
    avDiam(:,2) = mean(Y2,1,'omitnan');
    plot(permute(t, [2,1]), avDiam(:,2))
    hold off;
end

title('Average');
xlabel('Time (s)');
ylabel('Diameter (um)');
legend(RealID, ShamID);

baselineIdx = find(t(1,:) < stimOnset);
relAvDiam = avDiam./ mean(avDiam(1:max(baselineIdx),:),1,"omitnan")-1;
subplot(2,2,4);
plot(permute(t,[2,1]), relAvDiam(:,1)*100)
if isempty(Y2) == 0
    hold on
    plot(permute(t, [2,1]), relAvDiam(:,2)*100);
    hold off
end
title('Average - Relative');
xlabel('Time (s)');
ylabel('Diameter (%)');
legend(RealID, ShamID);


sgtitle(fig,VesselID);
savefig(fig,fullfile(out,strcat("diameterPlot_",VesselID)));

save(fullfile(out,strcat('diameters',VesselID,'.mat')),'bDiam','trialFlag','avDiam','relAvDiam','t','-mat','-append');