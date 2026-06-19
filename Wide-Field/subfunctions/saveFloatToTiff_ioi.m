function saveFloatToTiff_ioi(data,filename)
%% Description
% Saves a dF video as a 32-bit floating point tiff. 
%
% Written by Eric Martineau - Universtié de Montréal

%% Script
filename = strcat(filename, ".tif");
if isfile(filename) > 0
    delete(filename);
end

% Image must be single precision
if contains(class(data),'double') == 1
    data = single(data);
end

data(isnan(data)) = 0; %make NaN equal to 0 in case edges are NaN due to registration

% Create tiff object.
tiffObject = Tiff(filename, 'w');

% Set tags
tagstruct.ImageLength = size(data,1); 
tagstruct.ImageWidth = size(data,2);
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1; %for grayscale
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% Write the array to disk
tiffObject.setTag(tagstruct); 
tiffObject.write(data(:,:,1)); %write first image
for i = 2:size(data,3)  
    tiffObject.writeDirectory;
    tiffObject.setTag(tagstruct);
    tiffObject.write(data(:,:,i));    
end
tiffObject.close;
