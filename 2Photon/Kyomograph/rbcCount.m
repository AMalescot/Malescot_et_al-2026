function [flux, x, kmean, thresh,hThresh,fig] = rbcCount(kymo,kmean, thresh, hThresh, freq, dT_win, dT_step, X,discardLines)
%% Description:
% Script calculating RBC flux (rbc/s) from a kymograph by drawing a line
% across the center of the vessel and extracting the number of low points.

%% Input:
%  - kymo : stiched_linescan, must be T x X format (only 1 channel)
%  - path : self explanatory 
%  - freq : linerate of acquisition (in Hz)
%  - kmean : both the min peak distance for RBC count and the moving average window to smooth trace, in seconds
%    counted. Input empty vector if you want to determine manually
%  - thresh :  min peak prominence for RBC count. Input empty vector if you want to determine manually
%  - hThresh : max peak height for RBC count (max gray value to be considered a shadow). Input empty vector if you want to determine manually
%  - dT_win : delta-timewindow, in seconds. Determines the timewindow
%    over which the velocity is calculated
%  - dT_step : step by which the window is moved to calculate the velocity
%    over time. Must be equal or inferior to dT_win. Increases runtime if
%    too low.
%  - X : vector containing two x values for identifying vessel boundaries in
%    a line.
%  - discardLines : lost of line index to convert to exclude
 
%% Output:
%  - flux: RBC flux (in RBC/s)
%  - x: time axis (in s)
%  - kmean: both the min peak distance for RBC count and the moving average window to smooth trace, in seconds
%    counted. Will be used in the next trial
%  - thresh:  min peak prominence for RBC count. Will be used in the next
%    trial
%  - hThresh: max peak height for RBC count (max gray value to be considered a shadow). Input empty vector if you want to determine manually
%  - dT_win : delta-timewindow, in seconds. Will be used in the next
%    trial
%  - fig: figure 
%
% Script written by Éric Martineau - Université de Montréal
 %% Constants 
 dx = 3; % number of pixels each side of mid-vessel line to average to smooth signal
 
 %% Extract line across middle of vessel
mx = round(min(X)+(max(X) - min(X))/2);
line = kymo(:,(mx-dx):(mx+dx));

line = median(line,2);
 %% Smooth trace and query threshold 
 sLine = movmean(line,round(0.002*freq),'omitnan');
 sLine(discardLines) = NaN;


flag = 0;
if isempty(kmean) == 1
    flag = 1;
    answer = 'No';
    hThresh = mean(sLine,'omitnan');
    while contains(answer,'No')==1        
        fig = figure;
        x = 1000*((1:length(sLine))/freq);
        plot(x,sLine);
        xlabel('Time (ms)');
        ylabel('Pixel value');
        xlim([0 300]);
        grid on
        title('Min temporal difference between RBC shadows (in ms)');
        subtitle('Set low enough to distinguish two close RBC but high enough to not count one RBC as two');
        
        uiwait();
        prompt = {'Min Peak distance(ms)', 'Min Peak Prominence (Gray value)', 'Max peak height (Gray value)'};
        dlgtitle = "Processing parameters";
        dims = [1 50];
        definputs = {'2', '100', num2str(hThresh)};

        answer = inputdlg(prompt,dlgtitle, dims, definputs);

        kmean = str2double(answer{1})/1000; %threshold for rbc count
        thresh = str2double(answer{2});
        hThresh = str2double(answer{3});
        sLine = movmean(line,round(kmean*freq),'omitnan');
        sLine(discardLines) = NaN;

        %Validate
        [pk,RBCs] = findpeaks(-1*sLine,'MinPeakHeight',-hThresh, 'MinPeakDistance', kmean*freq,'MinPeakProminence',thresh);
        fig = figure('Position',[320,100,1120,840]);
        subplot(2,1,1);
        plot(sLine);
        xlabel('Time (lines)');
        ylabel('Pixel value');
        xlim([0 round(0.3*freq)]);%show 0.5seconds worth of lines
        grid on
        title('RBC counts');
        hold on 
        plot(RBCs,pk*-1, 'xg');
        hold off

        ax1 = subplot(2,1,2);
        axSize = round(ax1.Position(3:4).*fig.Position(3:4));
        im = imresize(kymo,[size(kymo,1) axSize(2)/(axSize(1)/(freq*0.3))]);
        imshow(mat2gray(rot90(im)),'Parent',ax1);
        xlim([0 round(0.3*freq)]);%show 0.3seconds worth of lines
        ax1.XAxis.Visible = 'on';
        title(strcat("Linerate = ", num2str(freq/10),"/100 ms"));
        xlabel('Time (lines)');
        uiwait();
        
        answer = questdlg('Is threshold OK?', ...
        'Validation', ...
        'Yes','No','No');
    end
end

%% Calculate bins
if flag == 0
    sLine = movmean(line,round(kmean*freq),'omitnan');
    sLine(discardLines) = NaN;
end
[pk,RBCs] = findpeaks(-1*sLine,'MinPeakHeight',-hThresh, 'MinPeakDistance', kmean*freq,'MinPeakProminence',thresh);

%Add a RBC in NaN intervals
indexes = diff(isnan(sLine)); %find changes to NaN values, 1 = number to NaN, -1 = NaN to number;
start = find(indexes == 1) + 1; %start index of nan intervals
finish = find(indexes == -1); %end index of nan intervals
middle = (finish + start)/2; %middle index of NaN interval
amps = mean(-1*[sLine(start-1) sLine(finish+1)],2); %average amplitude between NaN range, for graph purposes only
[RBCs,I] = sort([RBCs; middle]); %add indexes to RBC indexes and sort in ascending order
pk = [pk; amps]; %add amplitude and apply sort index from RBCs, for graphs only.
pk = pk(I);

lineWin = dT_win * freq; %number of lines in window
lineStep = dT_step * freq; %number of lines in step

flux = zeros(ceil(length(line)/lineStep),1);
x = zeros(ceil(length(line)/lineStep),1);

for i = 1:ceil(length(line)/lineStep)
    winStart = lineStep*(i-1)+1;
    winEnd = winStart + lineWin;
    winRBC = RBCs(RBCs > winStart);
    winRBC = winRBC(winEnd > winRBC);
    flux(i,1) = length(winRBC)/dT_win;
    x(i,1) = winEnd/freq;
end

%Rebuild invisible figure
fig = figure('Position',[320,100,1120,840]);
subplot(2,1,1);
plot(sLine);
xlabel('Time (lines)');
ylabel('Pixel value');
xlim([0 round(0.3*freq)]);%show 0.5seconds worth of lines
grid on
title('RBC counts');
hold on
plot(RBCs,pk*-1, 'xg');
hold off

ax1 = subplot(2,1,2);
axSize = round(ax1.Position(3:4).*fig.Position(3:4));
im = imresize(kymo,[size(kymo,1) axSize(2)/(axSize(1)/(freq*0.3))]);
imshow(mat2gray(rot90(im)),'Parent',ax1);
xlim([0 round(0.3*freq)]);%show 0.3seconds worth of lines
ax1.XAxis.Visible = 'on';
title(strcat("Linerate = ", num2str(freq/10),"/100 ms"));
xlabel('Time (lines)');
    
