function [y_fit, coefs,gof,out]= nl_detrend(x, data, bounds) 
% Fits a bi-exponential function, y = a*e^(x*b) + c*e^(x*b), on the average pixel, excluding the
% timepoints inside the bounds input (stim + expected reponse window) 
% Inputs:
%   x = time values as vector
%   data = video, inputed as XY-by-t matrix
%   bounds = exclusion bounds. Values to ignore during fit.
%   chanName = (optional) channel name, to display in waitbar. 
%
% Based on exemple from https://www.mathworks.com/matlabcentral/answers/431697-make-curve-fitting-faster
% Written by Éric Martineau and Antoine Malescot - Université de Montréal

%%
%Prep
if isempty(bounds)==0
    nPixels = size(data,1);
    binSize = 256;
    nBin = nPixels/binSize;
    nanIdx = isnan(data(1,:)); %identify nan frames (movement);
    data = data(:,~nanIdx); %remove NaN frames;
    xx = x(:,~nanIdx); %remove NaN timepoints, but keep original timescale for detrend function later
    indx = [find(xx>=bounds(1),1), find(xx>=bounds(2),1)];
    if length(indx)==2
        data = data(:,[1:indx(1),indx(2):end]); %remove points outside of bound
        xx_mat = repmat(xx(:,[1:indx(1),indx(2):end]),binSize,1);%remove points outside of bound
    elseif isempty(find(xx>=bounds(2),1))==1 %if no detrend points after stim, only detrend on before stim
        data = data(:,1:indx(1));
        xx_mat = repmat(xx(:,1:indx(1)),binSize,1);
    elseif isempty(find(xx>=bounds(1),1))==1 %if no detrend points before stim, only detrend on after stim
        data = data(:,indx(2):end);
        xx_mat = repmat(xx(:,indx(2):end),binSize,1);
    end

    %Fit average image to estimate starting parameters
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [1 -100 1 -0.1], 'Upper',[2^16 0 2^16 0.1],'TolX',1e-20,'TolFun',1e-20);
    opts.Display = 'Off';
    opts.StartPoint = [116.607174931279 -0.260820095436706 5036.74930225973 0.00129765068984856];
    [fitresult,gof,out] = fit( xx_mat(1,:)', mean(data,1)', 'exp2',opts);
    coefs = coeffvalues(fitresult); %coeff of mean pixel trend
    y_fit = fitresult(x);
    y_fit = y_fit';
else
    y_fit = ones(size(data,1),size(x,2));
    coefs = [];
    gof = struct();
    gof.sse = NaN;
    gof.rsquare = NaN;
    gof.adjrsquare = NaN;
    gof.rmse = NaN;
end

