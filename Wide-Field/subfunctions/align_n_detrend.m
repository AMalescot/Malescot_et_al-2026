function [data] = align_n_detrend(chan_name,infos,bounds)

%% DESCRIPTION %%
%Subfuction to simplify and streamline the preprocessing loop

%% Input:
% - chan_name: channel name.
% - infos: metadata for the LED channel.
% - bounds: exclusion bounds. Values to ignore during fit.
%
%% Output:
%  - data: detrended data
%
% Written by Éric Martineau and Antoine Malescot - Université de Montréal
%% Open a channel %%
fid = fopen(chan_name,'r+');
data = fread(fid, Inf, 'single=>double');
data = reshape(data,infos.datSize(1,1),infos.datSize(1,2),infos.datLength);  
data = reshape(data,[],infos.datLength);


%% Detrend channel %%
% Exclude stim and NaN frame and concatenate the rest
x = linspace(1/infos.Freq,infos.datLength/infos.Freq,infos.datLength); %timepoints

[T, ~,gof,out] = nl_detrend(x, data, bounds); %exponential fit on average video to extract trend in LED power changes. Can be done in pixelwise fashion but process time is much slower
T = repmat(T,[prod(infos.datSize) 1]); %repeat trend curve for every pixel as plotDetrend expects a pixelwise detrendingplotDetrend(x,data,T,chan_name,gof);
plotDetrend(x,data,T,chan_name,gof,out);
data = data./T;
data = reshape(data,infos.datSize(1,1),infos.datSize(1,2),infos.datLength);
clear T gof

%% Save data %%
data = reshape(data,[],1);
frewind(fid);
fwrite(fid, data ,'single');
fclose(fid);