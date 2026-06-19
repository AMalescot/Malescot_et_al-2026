function plotDetrend(x,data,T,chan_name,gof,out)
%% Description:
% function to plot the results of the detrending.
%% Input:
%  - x: array representing the timescale (1,data time length)
%  - data: array containing every pixel (256*256,data time length)
%  - T: trend calculated from bi-exponential function
%  - chan_name: name of the channel
%  - gof: quality metrics of the bi-exponential fitting
%
% Written by Éric Martineau - Université de Montréal
%% 
fig = figure;

%remove nanFrames for plot
nanFlag = isnan(data(1,:));

%Plot before detrending
subplot(2,1,1)
hold on
plot(x,mean(data,1),'-','Color',[0 0.4470 0.7410],'LineWidth',1)
plot(x,mean(T,1),'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1)

hold off
%Some useful info
title("Mean pixel value before detrending");
str = {strcat("SSE = ", num2str(gof.sse)), strcat("Rsquare = ", num2str(gof.rsquare)),strcat( "adjusted Rsquare = ", num2str(gof.adjrsquare)),strcat("RMSE = ", num2str(gof.rmse)),...
strcat("Iterations = ", num2str(out.iterations))};
str = [str strcat("Function eval = ", num2str(out.funcCount)) out.message];     

annotation('textbox', [0.55, 0.75, 0.1, 0.1], 'String', str,'FitBoxToText','on');

%Plot detrended data
subplot(2,1,2)
hold on
plot(x,mean(data./T,1),'-','Color',[0 0.4470 0.7410],'LineWidth',1);
hold off
title("Mean pixel value after detrending")

% Mark discarded areas
ax = fig.Children;
yy = [ax(1).YLim(1) ax(1).YLim(2) ax(1).YLim(2) ax(1).YLim(1)]; %y-limits of red zones
yy2 = [ax(2).YLim(1) ax(2).YLim(2) ax(2).YLim(2) ax(2).YLim(1)]; %y-limits of red zones
begin = 0; % 1 if begining of a series of nan values
for i = 1:length(x)
    if nanFlag(i) == 1 && begin == 0 %if encounter nan and not already in series
        xx = [x(i) x(i)]; %store begining of patch x-values
        begin = 1;
    elseif nanFlag(i) == 0 && begin == 1 %if encounter non-nan and in nan series
        xx = [xx x(i) x(i)]; %store end of patch x values
        patch(ax(1),xx,yy,'r','FaceAlpha',0.5,'EdgeColor','none'); %draw patch on first subplot
        patch(ax(2),xx,yy2,'r','FaceAlpha',0.5,'EdgeColor','none'); %draw patch on second subplot
        begin = 0; %mark exit of nan values series
    end %otherwise do nothing
end
sgtitle(chan_name);
savefig(fig,strcat("Detrend_",extractBefore(chan_name,'.dat'),".fig"),'compact');
close(fig)

