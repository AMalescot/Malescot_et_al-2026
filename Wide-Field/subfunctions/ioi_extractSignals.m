function ioi_extractSignals(AllSignals,Trials, ExtractID, SeriesID, PostStim, MaskThreshold, MaskChannel, avBin, root, Map, window, fps, Sig2Map)

%% Sig2Map = index of signals(stims)to map (1 = only first which is the real stims, 2 = both signals in case that two stim exist). All ROI signals will be extracted nevertheless
%% Filter out unwanted signals
AllSignals = AllSignals(ExtractID,:);
Trials = Trials(ExtractID);
SeriesID = SeriesID(ExtractID);

%% Load and log some parameters
chanIDX = fieldnames(AllSignals);
load(Trials{1,1});
StimStart = Trial(1).StimStart - Trial(1).Start;
StimDur = Trial(1).StimEnd - Trial(1).StimStart;
dimensions = size(AllSignals(1,1).Green);
FrameStart = round(fps*StimStart/1000);
FrameMaskEnd = round(fps*(StimStart+StimDur)/1000);
FrameEnd = round(fps*(StimStart+StimDur+PostStim)/1000);

%Map thresholds constants... could be changed in future uses
f1Threshold = 0;
f2Threshold = 0;
gThreshold = 0;
yThreshold = 0;
rThreshold = 0.5;

Parameters = struct('Masks', [],'MaskChannel', MaskChannel, 'MaskThreshold', MaskThreshold, 'gThreshold',gThreshold, 'yThreshold', yThreshold, 'rThreshold', rThreshold, 'f1Threshold', f1Threshold, 'f2Threshold', f2Threshold, 'RunAv', avBin,'RunAvFluo', avBin*2, 'PostStimAv', PostStim, 'FrameStart', FrameStart, 'FrameEnd', FrameEnd, 'StimStart', StimStart, 'StimDur', StimDur, 'windowMask', window, 'fps', fps, 'ExtractID', ExtractID, 'SeriesID',{SeriesID});

clear Trial StimStart StimDur frames ExtractID

%% Reorder signals %%
Spec_order = {"E2"; "D2"; "C2"; "B2"; "A2"; "E1"; "D1"; "C1"}; %%Specified order
idx = zeros(8,1);

for i = 1:length(idx)
     x = find(contains(SeriesID,Spec_order{i}));
     if isempty(x)==1
         idx(i,1) = NaN;
     else
         idx(i,1) = x;
     end
end
idx = rmmissing(idx);

AllSignals = AllSignals(idx,:);
SeriesID = SeriesID(idx,:);

clear idx Spec_order

%% Calculate Barrel masks from HbR signal from first stimulation %%
Masks = zeros(dimensions(1,1),dimensions(1,2),length(AllSignals));

for i = 1:size(AllSignals,1)
    if contains(MaskChannel,'F2')
        Signal = AllSignals(i,1).F2;
    elseif contains(MaskChannel,'F1')
        Signal = AllSignals(i,1).F1;
    elseif contains(MaskChannel,'Green')
        Signal = AllSignals(i,1).Green*-1;
    elseif contains(MaskChannel,'Amber')
        Signal = AllSignals(i,1).Amber*-1;   
    else
        Signal = AllSignals(i,1).Red*-1;
    end
    Signal = Signal .* window; %remove out of window signals before creating mask
    mask = BarrelMask(Signal,MaskThreshold, FrameStart, FrameMaskEnd,SeriesID(i,1)); %Filters barrel mask for area and eccentricity, see function code for parameter value
    Masks(:,:,i) = mask;
    clear HbR img mask i
end

Parameters.Masks = Masks; 

%% Map signals and extract ROI traces
mkdir('Analysis');
addpath('Analysis');
savepath
cd(strcat(root,'\Analysis'));

MaskMerge = zeros(size(Masks,1),size(Masks,2),3);

for i = 1:size(Masks,3)
    img = uint8(255 * mat2gray(Masks(:,:,i)));
    img = ind2rgb(img,custom_cmap(i));
    
    MaskMerge = max(MaskMerge, img);
    clear img
end

imwrite(MaskMerge,'Masks.jpeg');

for k = 1:size(AllSignals,2)
    filename = strcat("SignalMaps", num2str(k));
    save(filename, 'Parameters');
    
    for i = 1:5
        idx = chanIDX{i,1};
        switch idx
            case 'Green'
                A = [AllSignals(:,k).Green];
                if isempty(A) == 0
                    A = reshape(A, size(A,1),[],size(AllSignals,1), size(A,3));
                    A = -1*(A .* window);
                    if k <= Sig2Map
                        [gSig,gMap] = mapSignals(A, Map, gThreshold, FrameStart, FrameEnd);
                        save(filename, 'gSig', 'gMap', '-append');
                        imwrite(gSig,strcat('gSig', num2str(k),'.jpeg'));
                        imwrite(gMap,strcat('gMap', num2str(k),'.jpeg'));
                    end
                    [gtraces, fig] = pullROIsignal(A, Masks, avBin, fps, SeriesID, SeriesID);
                    sgtitle('HbT (green)');
                    save(filename, 'gtraces', '-append');
                    savefig(fig, strcat('gTraces',num2str(k)));                   
                                       
                end
                clear A fig

            case 'Amber'
                A = [AllSignals(:,k).Amber];
                if isempty(A) == 0
                    A = reshape(A, size(A,1),[],size(AllSignals,1), size(A,3));
                    A = -1*(A .* window);
                    if k <= Sig2Map
                        [ySig,yMap] = mapSignals(A, Map, yThreshold, FrameStart, FrameEnd);
                        save(filename, 'ySig', 'yMap', '-append');
                        imwrite(ySig,strcat('ySig', num2str(k),'.jpeg'));
                        imwrite(yMap,strcat('yMap', num2str(k),'.jpeg'));
                    end
                    [ytraces, fig] = pullROIsignal(A, Masks, avBin, fps, SeriesID, SeriesID);
                    sgtitle('HbT (yellow)');
                    save(filename, 'ytraces', '-append');
                    savefig(fig, strcat('yTraces',num2str(k)));     
                end
                clear A fig

            case 'Red'
                A = [AllSignals(:,k).Red];
                if isempty(A) == 0
                    A = reshape(A, size(A,1),[],size(AllSignals,1), size(A,3));
                    A = -1*(A .* window);
                    if k <= Sig2Map
                        [rSig,rMap] = mapSignals(A, Map, rThreshold, FrameStart, FrameEnd);
                        save(filename, 'rSig', 'rMap', '-append');
                        imwrite(rSig,strcat('rSig', num2str(k),'.jpeg'));
                        imwrite(rMap,strcat('rMap', num2str(k),'.jpeg'));
                    end
                    [rtraces, fig] = pullROIsignal(A, Masks, avBin, fps, SeriesID, SeriesID);
                    sgtitle('HbR (red)');
                    save(filename, 'rtraces', '-append');
                    savefig(fig, strcat('rTraces',num2str(k)));     
                end
                clear A fig

            case 'F1'
                A = [AllSignals(:,k).F1];
                if isempty(A) == 0
                    A = reshape(A, size(A,1),[],size(AllSignals,1), size(A,3));
                    A = A .* window;
                    if k <= Sig2Map
                        [f1Sig,f1Map] = mapSignals(A, Map, f1Threshold, FrameStart, FrameEnd);
                        save(filename, 'f1Sig', 'f1Map', '-append');
                        imwrite(f1Sig,strcat('f1Sig', num2str(k),'.jpeg'));
                        imwrite(f1Map,strcat('f1Map', num2str(k),'.jpeg'));
                    end
                    [f1traces, fig] = pullROIsignal(A, Masks, avBin*2, fps, SeriesID, SeriesID);
                    sgtitle('Fluo1 (475nm)');
                    save(filename, 'f1traces', '-append');
                    savefig(fig, strcat('f1Traces',num2str(k)));     
                end
                clear A fig

            case 'F2'
                A = [AllSignals(:,k).F2];
                if isempty(A) == 0
                    A = reshape(A, size(A,1),[],size(AllSignals,1), size(A,3));
                    A =A .* window;
                    if k <= Sig2Map
                        [f2Sig,f2Map] = mapSignals(A, Map, f2Threshold, FrameStart, FrameEnd);
                        [f2Sig2,f2Map2] = mapSignals(A, Map, 0.80, FrameStart, FrameEnd);
                        [f2Sig3,f2Map3] = mapSignals(A, Map, 0.75, FrameStart, FrameEnd);
                        [f2Sig4,f2Map4] = mapSignals(A, Map, 0.85, FrameStart, FrameEnd);
                        save(filename, 'f2Sig', 'f2Map', '-append');
                        imwrite(f2Sig,strcat('f2Sig', num2str(k),'.jpeg'));
                        imwrite(f2Map,strcat('f2Map', num2str(k),'.jpeg'));
                        imwrite(f2Sig2,strcat('f2Sig80_', num2str(k),'.jpeg'));
                        imwrite(f2Map2,strcat('f2Map80_', num2str(k),'.jpeg'));
                        
                        imwrite(f2Sig3,strcat('f2Sig75_', num2str(k),'.jpeg'));
                        imwrite(f2Map3,strcat('f2Map75_', num2str(k),'.jpeg'));
                        
                        imwrite(f2Sig4,strcat('f2Sig85_', num2str(k),'.jpeg'));
                        imwrite(f2Map4,strcat('f2Map85_', num2str(k),'.jpeg'));
                    end
                    [f2traces, fig] = pullROIsignal(A, Masks, avBin*2, fps, SeriesID, SeriesID);
                    sgtitle('Fluo2 (567nm)');
                    save(filename, 'f2traces', '-append');
                    savefig(fig, strcat('f2Traces',num2str(k)));     
                end
                clear A fig
        end
    end
end

clear idx gThreshold yThreshold rThreshold f1Threshold f2Threshold i k=