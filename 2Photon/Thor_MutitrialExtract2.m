%% Note: 
% Routine to align frames from multi trial acquisition.
% Data must be saved as a series of XY x T frames.

% The user needs to have saved one subfolder for each single trial
% containing metadata 'Experiment.xml' and data acquisition
% 'Image_scan_1_region_0_0.tif' saved in OME TIFF.
%
% On the main repository, a mat file has to be saved which is a structure
% contains all the individual trial information:
%             - Trial.Type: string identifier of the trial type ("real" or
%             "sham").
%             - Trial.Start: time of the trial start in us.
%             - Trial.StimStart: time of the stimulation start, in us.
%             - Trial.StimEnd: time of the stimulation end, in us.
%             - Trial.End: time of the end of the trial, in us.

%  Written by Eric Martineau and Antoine Malescot - Universite de Montreal
tic
%% Extract parameters
alignChan = 2; % Channel used to calculate alignment vectors
discardTrials = []; % Trial to discard
tBin = 3; %Number of frames to bin across time to increase SNR.
oemTiff = 1; % Value of 0 if images saved in raw tiff (series of images), value of 1 if saved in OEM tiff format (one image contains stack)
Stim = 1; % Value of 0 if there is no stim (no trial file), value of 1 if there is a trial(stack file)
align = 1; % Value of 1 if you want to align frames  within stack(t-stack), value of 0 if you don't want to (z-stack).

%% Import image files and trial data
clear RealData ShamData avRealData avShamData
ok = 0;
longAcq = 0;
while ok == 0
    % Get file path
    path = uigetdir(pwd,'Select image folder');
    parts = regexp(path,'\','split');

    %%% Check Series format %%%
    cd(path) % Open path to experiment repository
    D = dir; % Open struct with what it contains
    D = D([D.isdir]); % Check for subfolders  
    D = D(~contains({D.name},'.'));
    D = D(~contains({D.name},'jpeg')); % Keep the one with Series inside
    D = D(~contains({D.name},'AQuA'));
    D = D(~contains({D.name},'Vasometrics'));
    D = D(~contains({D.name},'Neuronal_activity'));
    
    if isempty(D)==1
        multitrial = 0;
        type = 'Single trial';
        nTrial = 1; 
        cd(fullfile(parts{1:end-1}));
        D = dir;
        D = D(contains({D.name},parts{end}));
        cd(path)
    else
        multitrial = 1; 
        type = 'Multi-trial';
        nTrial = length(D);        
    end

    %%% Check with users %%%
    prompt = {[strcat('Analyzing "', parts{end},'" ') newline strcat(type,' acquisition.')]};
    answer = questdlg(prompt,'Folder check','Yes','No - Reselect','Yes');

    if contains(answer,'Yes')==1
        ok = 1;
    else
        ok = 0;
    end
end

%%% Load trial file %%%
if Stim == 1
    [trialFile, trialPath] = uigetfile('.mat','Select Trial file');
    load(fullfile(trialPath, trialFile),'Trial');

    if nTrial == 1 && multitrial == 1
        nStim = length(Trial);
        longAcq = 1;
    end
end

%%% Remove discarded trials %%%
D(discardTrials) = [];
if Stim == 1
    Trial(discardTrials) = [];
    nTrial = length(D);
end


%%% Create Data file structure %%%
AllData = struct('ID',cell(nTrial,1),'Frames', cell(nTrial,1),...
    'Type',cell(nTrial,1),'Rewarded',cell(nTrial,1),...
    'metadata',cell(nTrial,1),'xmlMeta',cell(nTrial,1));
M = cell(1,nTrial);

%%% Import and alignement loop %%%
parfor k = 1:nTrial
    % Import
    if oemTiff == 1
        [stack, metadata, xmlMeta] = ThorImport(fullfile(D(k).folder,D(k).name), 'Image_scan_1_region_0_0.tif');
    elseif oemTiff == 0
        [stack, metadata, xmlMeta] = ThorImport_noOEM(fullfile(D(k).folder,D(k).name));
    end

    % Temporal binning
    if tBin > 1
        stack = TempBin(stack, tBin);
    end
    metadata.tBin = tBin;
    metadata.FrameRate = metadata.FrameRate/tBin;

    % Calculate within trial movement (rigid)
    if align == 1
        disp(strcat("Aligning stack#", num2str(k)," ...."));
        f_id = (1:size(stack,3)); % Index of each frame used by the alignment algorithm
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.GradientMagnitudeTolerance = 1.0e-04; 
        optimizer.MaximumStepLength = 0.00625;
        optimizer.MaximumIterations = 100; 
        optimizer.RelaxationFactor = .5; 
        fixedRefObj = imref2d(size(stack(:,:,1,alignChan)));

        [displace,m] = recursiveAlign_rigid(stack(:,:,:,alignChan),f_id,optimizer,metric, fixedRefObj); % Calculation of displacement vectors for each frame.

        %Apply transform
        for i = 1:size(stack,3) % Applies displacement field to each frame
            stack(:,:,i,:) = imwarp(stack(:,:,i,:),displace{1,i},'Outputview',fixedRefObj);
            disp(strcat("Aligning frame ",num2str(i),"/",num2str(size(stack,3))));
        end
        M(1,k) = {m};
    else
        disp(strcat("within trial alignement skipped for trial", num2str(k), "...."));
    end

    %Store data
    if longAcq == 0
        if multitrial == 1
            AllData(k).ID = D(k).name;
        else
            AllData(k).ID = parts{end};
        end
        AllData(k).Frames = stack;
        AllData(k).metadata = metadata;
        AllData(k).xmlMeta = xmlMeta;

        % Import relevant trial data
        if Stim == 1
            AllData(k).Type = Trial(k).Type;
            AllData(k).Rewarded = Trial(k).Rewarded;
        else
            AllData(k).Type = "Real";
        end
    else
        AllData(k).ID = "Long";
        AllData(k).Frames = stack;
        AllData(k).metadata = metadata;
        AllData(k).xmlMeta = xmlMeta;
    end
end
clear k

if longAcq == 1
    startOffset = Trial(1).Start;
    trialLength = ceil((Trial(1).End - startOffset)/1000 * AllData(1).metadata.FrameRate);
    startIdx = zeros(nStim,1);
    endIdx = zeros(nStim,1);
    stimIdx = zeros(nStim,2);

    for i = 1:nStim
        Trial(i).Start = Trial(i).Start - startOffset; % Start time relative to acquisition start
        Trial(i).StimStart = Trial(i).StimStart - startOffset;% Stim start time relative to acquisition start
        Trial(i).StimEnd = Trial(i).StimEnd - startOffset;% Stim end time relative to acquisition start
        startIdx(i) = ceil(Trial(i).Start/1000 * AllData(1).metadata.FrameRate)+1;
        endIdx(i) = startIdx(i) + trialLength-1;
        stimIdx(i,:) = [ceil(Trial(i).StimStart/1000  * AllData(1).metadata.FrameRate)+1, ceil(Trial(i).StimEnd/1000  * AllData(1).metadata.FrameRate)];

        % Store acquisition
        AllData(i+1).ID = num2str(i);
        AllData(i+1).Frames = AllData(1).Frames(:,:,startIdx(i):endIdx(i),:);
        AllData(i+1).metadata = AllData(1).metadata;
        AllData(i+1).xmlMeta = AllData(1).xmlMeta;
        AllData(i+1).Type = Trial(i).Type;
        AllData(i+1).Rewarded = Trial(i).Rewarded;
    end
end

% Perform between-trial alignments
if nTrial > 1
    if align == 0
        for k = 1:nTrial
            M{1,k} = max(AllData(k).Frames(:,:,:,alignChan),[],3);
        end
    end
    m1 = M{1,1}; % Fixed average trial1 image

    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.00625;
    fixedRefObj = imref2d(size(m1));

    parfor k = 2:nTrial
        tform = imregtform(M{1,k},m1,'rigid', optimizer, metric); % Calculate transform
        disp(strcat("Between trial alignement for stack#", num2str(k)," ...."));
        AllData(k).Frames = imwarp(AllData(k).Frames,tform,'Outputview',fixedRefObj);
    end
end

% Save everything
if Stim == 1
    fname = strcat(extractBefore(trialFile,'stack'),'AllData');
else
    fname = strcat(parts{end},'_AllData');
end
save(fullfile(path,fname),'AllData','discardTrials','alignChan','-v7.3');

%% Average trials %%
if multitrial == 1
    if longAcq == 1
        AllData(1) = [];
    end
    NbChan = AllData(1).metadata.NbChannels;
    flag = contains([AllData.Type],'Real');
    RealTrial = AllData(flag);
    ShamTrial = AllData(~flag);

    ShamData = [];
    RealData = [];

    for i = 1:length(RealTrial)
        RealData = cat(5,RealData,RealTrial(i).Frames);
    end

    avRealData = mean(RealData,5);
    for j = 1:NbChan
        saveTifStack(uint16(avRealData(:,:,:,j)),fullfile(path,strcat('avRealData_Ch',num2str(j))));
    end
    clear RealTrial

    if isempty(ShamTrial) == 0
        for i = 1:length(ShamTrial)
            ShamData = cat(5,ShamData,ShamTrial(i).Frames);
        end

        avShamData = mean(ShamData,5);
        for j = 1:NbChan
            saveTifStack(uint16(avShamData(:,:,:,j)),fullfile(path,strcat('avShamData_Ch',num2str(j))));
        end
    end
    clear ShamTrial
end
toc