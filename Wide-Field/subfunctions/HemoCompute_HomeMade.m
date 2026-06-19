function [HbO, HbR] = HemoCompute_HomeMade(DataFolder, SaveFolder, Parameters, Illumination, b_normalize)
%% Description:
% function approximating concentration variation of oxygenated (HbO) and
% de-oxygenated (HbR) hemoglobin from two or three illumination wavelengths
% of intrinsic signals.
%
%% Input:
%   - DataFolder (char): path to folder where the input files "red",
%     "green" are stored.
%   - SaveFolder (char): path where to save the HbO and HbR files. If empty,
%     the data will not be saved to a .dat file.
%   - FilterSet (char): type of excitation/emission set of filters used. It is one
%     of the following:
%           - gCaMP
%           - jrGECO
%           - none
%   - Illumination (cell): names of illumination colors used ("red", "green"). A minimal of two must be provided.
%   - b_normalize (bool): Set to TRUE, to normalize the data (when the input
%     data is not normalized already).
%% Output:
%   - HbO (numerical array): data with the same dimensions of the reflectance
%     channels containing the approximate variations of the oxygenated
%     hemoglobin.
%   - HbR (numerical array): data with the same dimensions of the reflectance
%     channels containing the approximate variations of the deoxygenated
%     hemoglobin.
%
% Original script from the UMIT toolbox. 
% Modified by Eric Martineau and Antoine Malescot to add other
% filter-wavelength combinations.

%% Preparation %%

FilterSet = Parameters.FilterSet;
Calcium = Parameters.Calcium;

%Inputs Validation
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = strcat(DataFolder, filesep);
end
bSave = true;
if isempty(SaveFolder)
    bSave = false; % Disable saving if SaveFolder was not provided.
elseif( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = strcat(SaveFolder, filesep);
end

if( ~contains(lower(FilterSet), {'vglut2_chr2_jrgeco1a','gcamp', 'jrgeco', 'none'}) )
    disp('Invalid Filter set name');
    return;
end

if( sum(contains({'red', 'green'}, lower(Illumination))) < 2 )
    disp('At least two different illumination wavelengths are needed for Hb computation');
    return;
end

% Tags:
fTags = {'fidR', 'fidG'};
cTags = {'iRed', 'iGreen'};
colors = {'Red', 'Green'};

%Files Opening:
NbFrames = inf;
fidR = 0;
fidG = 0;
for indC = 1:size(Illumination,2)
    switch lower(Illumination{indC})
        case 'red'
            fidR = fopen([DataFolder 'red.dat']);
            iRed = matfile([DataFolder 'red.mat']);
            NbFrames = min([NbFrames, iRed.datLength(1,end)]);
            NbPix = iRed.datSize;
            Freq = iRed.Freq;
        case 'green'
            fidG = fopen([DataFolder 'green.dat']);
            iGreen = matfile([DataFolder 'green.mat']);
            NbFrames = min([NbFrames, iGreen.datLength(1,end)]);
            NbPix = iGreen.datSize;
            Freq = iGreen.Freq;

        otherwise
            disp('Unknown colour');
    end
end

% Get data size:
for i = 1:2
    if exist(cTags{i},'var')
        break
    end
end
eval(['iFile = ' cTags{i} ';']);
datsz = [iFile.datSize, iFile.datLength];
datsz(end) = NbFrames; % Update trial/movie length.
NbPix = double(NbPix);


%% Check if the data is normalized: %%
indxNorm = [-2 -2];
disp('Checking channel data...')
for i = 1:2
    if eval([fTags{i} '== 0'])        
        continue
    end
    eval(['Mdat = mean(fread(' fTags{i} ', Inf, ''*single''),''all'', ''omitnan'');']);
    if Mdat >.75 && Mdat <1.25
        % If the data is centered at one.
        indxNorm(i) = 1;
    elseif Mdat > -.25 && Mdat <.25
        % If the data is centered at zero
        indxNorm(i) = 0;
    else
        % Data not normalized!
        indxNorm(i) = -1;
    end
end

%% Check if all channels have the same profile: %%
indxNorm(indxNorm == -2) = []; %remove -2 values
if ~all(indxNorm == indxNorm(1)) %homemade correction
    error('The input data is heterogeneous! All channels must be preprocessed in the same way.')
end
indxNorm = indxNorm(1);
if indxNorm == -1 && ~b_normalize
    error('Operation aborted! The channels must be normalized or set "b_normalize" input to TRUE.')
end
if indxNorm == 0
    warning('The channels are centered at zero. They will be shifted to be centered at one.')
end
if b_normalize && (indxNorm == 1 || indxNorm == 0)
    b_normalize = false;
    warning('The input data is already normalized. Normalization will be skipped!')
end
disp('Data checked!')

%% Filter setting %%
switch( lower(FilterSet) )
    case 'vglut2_chr2_jrgeco1a'
        Filters.Excitation = 'VGLut2_ChR2_jRGECO1a';
        Filters.Emission = 'VGLut2_ChR2_jRGECO1a';
    case 'gcamp'
        Filters.Excitation = 'GCaMP';
        Filters.Emission = 'GCaMP';
    case 'jrgeco'
        Filters.Excitation = 'jRGECO';
        Filters.Emission = 'jRGECO';
    otherwise
        Filters.Excitation = 'none';
        Filters.Emission = 'none';
end
Infos = load([DataFolder 'AcqInfos.mat']);
Filters.Camera = Infos.AcqInfoStream.Camera_Model;
clear Infos;

%% Computation itself %%
A = ioi_epsilon_pathlength_RL('Hillman', 100, 60, 40, Filters,Calcium);
f = fdesign.lowpass('N,F3dB', 4, 1, Freq); %Low Pass
lpass_high = design(f,'butter');
f = fdesign.lowpass('N,F3dB', 4, 1/120, Freq); %Low Pass
lpass_low = design(f,'butter');
if numel(datsz) == 4 % For 4D data with dimensions {'E','Y','X','T}:
    
    nIter = double(datsz(1));
    offset = 4;
    Size = [prod(datsz(2:3)), datsz(4)];
    Precision = '*single';
    Skip = (datsz(1)-1)*4;
else % For 3D data with dimensions {'Y','X','T}:
    MemFact = 16;
    Precision = [int2str(NbPix(1)*MemFact) '*single'];
    Skip = (NbPix(1)*NbPix(2) - NbPix(1)*MemFact)*4;
    nIter = NbPix(2)/MemFact;
    offset = NbPix(1)*MemFact*4;
    Size = [NbPix(1)*MemFact, NbFrames];
end

HbO = zeros(datsz, 'single');
HbR = zeros(datsz, 'single');

%% Computation loop %%
h = waitbar(0,'Computing');
for indP = 1:nIter
    if( fidR )
        fseek(fidR,(indP-1)*offset,'bof');
        Red = fread(fidR,Size,Precision,Skip);
        if b_normalize
            Red = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Red)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Red)'))';
            tmp(tmp<min(Red(:))) = min(Red(:));
            Red = (Red)./(tmp);
        end
        if indxNorm == 0 % data centered at zero
            Red = Red + 1;
        end
        Red = -log(Red);
    end
    if( fidG )
        fseek(fidG,(indP-1)*offset,'bof');
        Green = fread(fidG,Size,Precision,Skip);
        if b_normalize
            Green = single(filtfilt(lpass_high.sosMatrix, lpass_high.ScaleValues, double(Green)'))';
            tmp = single(filtfilt(lpass_low.sosMatrix, lpass_low.ScaleValues, double(Green)'))';
            tmp(tmp<min(Green(:))) = min(Green(:));
            Green = (Green)./(tmp);
        end
        if indxNorm == 0 % data centered at zero
            Green = Green + 1;
        end
        Green = -log(Green);
    end
    clear tmp;

    Ainv = pinv(A(1:2,:));
    Hbs = Ainv*([Red(:), Green(:)]') .* 1e6;
    clear Red Green;


    if numel(datsz) == 4
        Hbs = reshape(Hbs, [2 1, datsz(2:end)]);
        Hbs = real(Hbs);
        HbO(indP,:,:,:) = squeeze(Hbs(1,:,:,:,:));
        HbR(indP,:,:,:) = squeeze(Hbs(2,:,:,:,:));
    else
        Hbs = reshape(Hbs, 2, NbPix(1), MemFact, []);
        Hbs = real(Hbs);
        HbO(:,(indP-1)*MemFact + (1:MemFact),:) = squeeze(Hbs(1,:,:,:));
        HbR(:,(indP-1)*MemFact + (1:MemFact),:) = squeeze(Hbs(2,:,:,:));
    end

    waitbar(indP/nIter,h);
end
close(h);

%% Save File management %%
if( bSave )
    % Save HbO:
    % Save .DAT file:
    fidHbO = fopen([SaveFolder 'HbO.dat'],'W');
    fwrite(fidHbO, HbO, '*single');
    fclose(fidHbO);
    fn = setdiff(fieldnames(iFile), {'Properties', 'datFile'});

    % Save .MAT file:
    fHbO = matfile([SaveFolder 'HbO.mat'], 'Writable', true);
    fHbO.datFile = 'HbO.dat';
    for i = 1:numel(fn)
        fHbO.(fn{i}) = iFile.(fn{i});
    end

    % Save HbR:
    % Save .DAT file:
    fidHbR = fopen([SaveFolder 'HbR.dat'],'W');
    fwrite(fidHbR, HbR, '*single');
    fclose(fidHbR);

    % Save .MAT file:
    fHbR = matfile([SaveFolder 'HbR.mat'], 'Writable', true);
    fHbR.datFile = 'HbR.dat';
    for i = 1:numel(fn)
        fHbR.(fn{i}) = iFile.(fn{i});
    end
end

end