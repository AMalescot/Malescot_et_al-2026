clear all
close all
clc
%% Note:
% Here is a step by step preprocessing for wide field optical imaging data.
% This program requires several functions from the UMIT github (https://github.com/LabeoTech/Umit).
%% Initialization
root = 'G:\sandbox_WF';

%%% Parameters %%%
Parameters.HemoCorrChan = {'Green&Red'};
Parameters.Sigma = 3; % Value for gaussian filter
Parameters.base_perc = 0.75; % Baseline proportion to keep, value between 0 and 1
Parameters.tBin = 1; % Temporal binning value
Parameters.sBin = 1; % Spatial binning value
Parameters.final_tBin = 4; % Final framerate for video preview
Parameters.ComputeHbO = 1; % HbO computation (boolean value)
Parameters.alpha = 1; % Weight of the excitatory component in beer-lambert law
Parameters.beta = 1; % Weight of the emission component in beer-lambert law
Parameters.FilterSet = 'jRGECO'; % Char vector defining filterset ('jRGECO','gcamp','vglut2_chr2_jrgeco1a','none'). Refers to filter configurations defined in ioi_epsilon_pathlength_RL(). The filter names and spectrums are loaded from SysSpect_RL.mat. 
Parameters.Calcium = 'jRGECO1a'; % Calcium indicator
Parameters.DetrendBounds = [4, 10]; % Time interval to exclude from detrending
Parameters.Discard = [];

hemoCorr = {'Green','Red'};
saveRaw = 0; % Save raw .dat and data.mat files (boolean value)
indivTrial = 0; % Save single trial data (boolean value)

%%% Extract other useful path %%%
root_serie = fullfile(root,'Series1'); % Path to series experiment
root_trial = fullfile(root,'19-Oct-2023_#3501_VGLUT2_ChR2_jRGECO1a_stack1.mat'); % Path to series experiment
root_img = fullfile(root, 'Snapshot20231019_124159_Green.png'); % Path to snapshot

load(root_trial) % load Trial file

%% Step 1: Pre-processing %%
[SignalsFilt, Parameters.fps] = ioi_preprocessing(root_serie, Trial, Parameters, hemoCorr, saveRaw, indivTrial);
SignalsFilt(1).HbT = SignalsFilt(1).HbO + SignalsFilt(1).HbR;
SignalsFilt(2).HbT = SignalsFilt(2).HbO + SignalsFilt(2).HbR;

save(fullfile(root,'Processed_Signals.mat'),'SignalsFilt','-v7.3')
mkdir(fullfile(root,'Analysis'))
save(fullfile(root,'Analysis','Parameters.mat'),'Parameters')

%% Step 2: Find activation mask %%

%%% Parameters %%%
Parameters.MaskChannel = 'F2'; % Reference channel to extract mask
Parameters.MaskThreshold = 0.85; % Percentage to keep in normalized signal
Parameters.FrameStart = 80; % First frame of the average
Parameters.FrameEnd = 180; % Last frame of the average

Parameters.windowMap = imread(root_img);

%%% Draw Window Mask %%%
figure()
hold on
title('Draw window mask')
imshow(Parameters.windowMap)
roi = images.roi.AssistedFreehand;
draw(roi);
Parameters.windowMask = createMask(roi);
imshow(Parameters.windowMap.*uint8(Parameters.windowMask))
close all

%%% Extract ROI %%%
[Parameters.Mask] = Signal_Mask(SignalsFilt(1).(Parameters.MaskChannel), Parameters.MaskThreshold, Parameters.FrameStart, Parameters.FrameEnd);
save(fullfile(root,'Analysis','Parameters.mat'),'Parameters')

%% Step 3: Extract Signal %%

%%% Parameters %%%
Parameters.RunAv = 5; % Moving average
[traces] = pullROIsignal(SignalsFilt(1), Parameters.Mask, Parameters.RunAv);
save(fullfile(root,'Analysis','SignalMaps.mat'),'traces')
save(fullfile(root,'Analysis','Parameters.mat'),'Parameters')

%% Step 4: Donut analysis %%

%%% Parameters and options %%%
options.iterations = 10; % Number of donuts
options.radius = 5;  % Initial mask radius
[Signal,Mask_donut] = Analysis_donuts_timbitsApp(SignalsFilt,Parameters,options);
save(fullfile(root,'Analysis','Donut_Analysis.mat'),'Signal','Mask_donut')