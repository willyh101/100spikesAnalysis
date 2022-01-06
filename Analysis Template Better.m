clear;
date = '200727';
mouse = 'W26_1';%'I138_1';%'I136_1';
epochs = '2_3_4_5_7_8_9_10_12';

% addpath(genpath('C:\Users\Will\Lab Code\Ian Code'))
% basePath = ['E:\Contrast Modulated Ensembles\' mouse '\' date '\'];
basePath = ['C:\Users\ian\Documents\DATA\F\' mouse '\' date '\'];
% basePath = ['T:\Ian\F\' mouse '\' date '\'];

path = fullfile(basePath,epochs);

baseName = [mouse '_' date];%'I118a.2_180504';
loadList = {['F_' baseName '_plane1_proc'] ['F_' baseName '_plane2_proc'] ['F_' baseName '_plane3_proc']};% ['F_' baseName '_plane4_proc']};

nDepthsTotal = 3;4;%Normally 3;
physfile = fullfile(basePath,[date '_A' '.mat']);
% physfile = fullfile(basePath,[date(3:end) '_A' '.mat']);
disp('Loading...')
try
    load(physfile)
catch
    physfile = fullfile(basePath,[date(3:end) '_A' '.mat']);
    load(physfile)
end
disp('Loaded')
%% Experiment
%
% s2pEpoch = 2 ;
% DAQepoch = 2 ;
theList=[];

%order s2pEpoch DAQepoch out condition 
% condition options are 'stim' 'exp' 'vis' 'vis2' 'exp2' 'mani' 'spk' or
% 'info' ('info' is included in 'exp' but can also be overwritten alone)
% spk is an extra module on exp, so run exp first even if it will be
% overrun
theList = {
    3 4 'stim'
    5 7 'vis2'
    7 9 'vis'
    8 10 'exp'
    8 10 'spk'
    9 12 'exp'
    9 12 'mani'
    ...6 6 'exp2'
    ...5 5 'info'
    };
listSize = size(theList);



%% Scary Loading Part
%The Slow File reading Part

for listEntry = 7:listSize(1);
    s2pEpoch = theList{listEntry,1};
    DAQepoch =  theList{listEntry,2};
    option =  theList{listEntry,3};
    
    disp(['Processing ' option '. s2p: ' num2str(s2pEpoch) '. DAQ: ' num2str(DAQepoch)]);
    
    BaseLinePeriod = 1000;
    
    if ~exist('ExpStruct')
        disp('Reloading ExpStruct')
        load(physfile,'ExpStruct');
    end
    
    if ~exist('DIGITSWEEP')
        DIGITSWEEP = ExpStruct.digitalSweeps;
    end
    
    numColors=2;
    
    tScary = tic;
    [localFiles, frameCount, MD, FR, minNumFrames, numFrames, numCells,...
        numFiles, allData, allSPData, allCoM, allDepth, allNP, allNPC,...
        roiMasks, T, BL, allDataCorrect, zData, zSPData, backupData,...
        DFFdata, backupSPData, ExpStruct, CellstToExclude, allDataNoNP,...
        stimID, uniqueStims, uStimCount, runVector, ArtifactSizeLeft,...
        ArtifactSizeRight, ROIinArtifact, excludeCells, framesUsed, Tused, zData2]...
        = loadCaData(path,loadList,DAQepoch,s2pEpoch,nDepthsTotal,ExpStruct,DIGITSWEEP,BaseLinePeriod, numColors);
    toc(tScary)
    
    swp =ExpStruct.EpochEnterSweep{DAQepoch};
    hr = ExpStruct.Holo.Sweeps_holoRequestNumber(swp);
    
    holoRequests = ExpStruct.Holo.holoRequests{hr};
    
    if ~isfield(holoRequests,'rois')
        newHR = ExpStruct.Holo.holoRequests{hr-1}; %catch error where sometimes wrong
        
        if isfield(newHR,'rois')
            disp('HR.rois not detected; replacing HR with newHR')
            holoRequests = newHR;
        end
    end
    
    
    [dfData, zdfData] =  computeDFFwithMovingBaseline(allData);
    
    
    % %%
    % backupStimID = stimID;
    % backupRunVector = runVector;
    %%Experiment Design
    
    strt = BaseLinePeriod; %start Time;
    try
        numStimuli = holoRequests.holoStimParams.nHolos(1); %number of stimululations per trial. (c
        freqOfStim = holoRequests.holoStimParams.hzList(1)/holoRequests.holoStimParams.pulseList(1); %repetition rate of stimuli in Hz.
    catch
        disp('***holoStimParams Error***')
        numStimuli=1;
        freqOfStim=1;
    end
    stimulusLength = 1/freqOfStim*1000; %how long stimulating for in ms
    
    
    numHolosPerTrial = 1; %numStimuli;
    dura = stimulusLength;
    freq = freqOfStim;
    
    %%Determine runTrials, Motion Trials
    if ~exist('ExpStruct')
        load(physfile,'ExpStruct');
    end
    
    if ~exist('numFiles')
        numFiles = numel(frameCount);
    end
    
    if ~exist('runVector')
        [runVector] = runReader2(ExpStruct,DAQepoch,FR,3);% sometimes 3, sometimes 5
    end
    
    numStimuli = 1;
    
    strtFrame = round(strt/FR);
    
    motionThreshold = 3;5;  %how many pixels can you offset before trial excluded -important for holography
    runThreshold = 10; % 6 or 10 standard; Run Speed requirement (evan says 6, i was doing 10 before).
    
    
    % lowMotionTrials = false([3 numTrials numFiles]);
    lowMotionTrials = false([numStimuli numFiles]);
    highRunTrials = false([numStimuli numFiles]);
    
    Toffsets = cellfun(@(x) median(x),Tused,'uniformoutput',0);
    
    tTest = cellfun(@(x) x-median(x),Tused,'uniformoutput',0);
    
    Tarray = cellfun(@(x) reshape(x,numFrames,numFiles,2),tTest,'uniformoutput',0);
    runVal=[];
    for i = 1:numStimuli
        
        idx = max(round((strt/1000+(i-1)/freqOfStim)*FR),1):min(round((strt/1000+(i-1)/freqOfStim)*FR)+stimulusLength/1000*FR,size(runVector,2));
        for k = 1:numFiles
            m = cellfun(@(x) squeeze(mean(x(idx,k,:),1)), Tarray,'uniformoutput',0);
            m = abs(cell2mat(m));
            m = mean(m,2);
            mLow = all(m<motionThreshold);
            lowMotionTrials(i,k)=mLow;
            highRunTrials(i,k)=mean(runVector(k,idx)>runThreshold)>0.5;
            runVal(i,k)= mean(runVector(k,idx));
        end
    end
    
    figure(1002);
    subplot(1,3,1)
    imagesc(lowMotionTrials)
    title('Low Motion Trials')
    subplot(1,3,2)
    imagesc(highRunTrials)
    title('High Run Trials')
    subplot(1,3,3)
    imagesc(runVector')
    %%Align Targets Shot to Cells
    
    
    numStims = size(holoRequests.targets,1);
    depthList = unique(holoRequests.targets(:,3));
    
    stimCoM=[];
    stimDepth=[];
    
    offsets = -mean(cat(1,Toffsets{:}));
    if ~exist('dat')
        load(fullfile(path,[loadList{1} '.mat']),'dat');
    end
    
    additionalOffsets = [0 0 ];
    offsets = offsets + additionalOffsets;
    
    
    if isfield(dat.ops,'cutArtifact') && dat.ops.cutArtifact
        offsets = offsets+ [0 -dat.ops.ArtifactStart];
    end
    % stimCoM = fliplr(holoRequests.targets(:,1:2))+offsets;
    stimCoM = holoRequests.targets(:,1:2)+offsets;
    stimDepth = holoRequests.targets(:,3);
    
    for i = 1:numel(depthList)
        stimDepth(stimDepth==depthList(i))=i;
    end
    
    Mapping =zeros([numStims 3]);
    for s = 1:numStims
        
        d = stimDepth(s);
        Mapping(s,1) = d;
        
        a = allCoM;
        a(allDepth~=d,:) = nan;
        
        b = stimCoM(s,:);
        [Mapping(s,3), Mapping(s,2)] = min ( sqrt(sum((a-b).^2,2)) );
    end
    
    targetDistanceThreshold = 10; %has been 15
    targettedCells = Mapping(:,2); targettedCells(Mapping(:,3)>targetDistanceThreshold)=nan;
    targettedCells(ismember(targettedCells,CellstToExclude))=nan;
    
    figure(1);clf
    for i = 1:nDepthsTotal
        subplot(1,nDepthsTotal,i)
        plot(allCoM(allDepth==i,1),allCoM(allDepth==i,2),'o')
        hold on
        plot(stimCoM(stimDepth==i,1),stimCoM(stimDepth==i,2),'.')
        xlim([0 512])
        ylim([0 512])
    end
    
    % % new method using KDTree
    % % this can be done in 3D, so make the depths correct
    % allDepthsNew = arrayfun(@(x) depthList(x), allDepth);
    % a = [allCoM, allDepthsNew];
    % b = [stimCoM, holoRequests.targets(:,3)];
    %
    % % make kdtree and query
    % kmodel = KDTreeSearcher(b);
    % [targettedCells, dists] = knnsearch(kmodel, a);
    % targettedCells(dists > targetDistanceThreshold)=nan;
    % targettedCells(ismember(targettedCells,CellstToExclude))=nan;
    
    
    % targettedCells = Mapping(:,2); targettedCells(Mapping(:,3)>targetDistanceThreshold)=nan;
    % targettedCells(ismember(targettedCells,CellstToExclude))=nan;
    disp(['Number of missed targets: ' num2str(sum(isnan(targettedCells)))]);
    disp(['Thats ' num2str(sum(isnan(targettedCells))/numel(targettedCells)*100) '% of Targets were missed'])
    
    
    %%New ROIinArtifact
    
    ArtifactSizeLeft = 100;
    ArtifactSizeRight = 100;
    
    yoffset = -offsets(2);
    % ROIinArtifactBackup = ROIinArtifact;
    ROIinArtifact = allCoM(:,2)<ArtifactSizeLeft-yoffset | allCoM(:,2)>511-(ArtifactSizeRight+yoffset);
    
    %%Determine holograms shot
    try
        roisTargets = holoRequests.rois;
        
        HoloTargets=[];
        
        numHolos = numel(roisTargets);
        
        for i = 1:numHolos
            HoloTargets{i} = targettedCells(roisTargets{i})';
        end
        
        try
            TargetRois = unique(cat(2,roisTargets{:}));
        catch
            TargetRois = unique(cat(1,roisTargets{:}));
        end
        allRoiWeights = holoRequests.roiWeights;
        TargetRoiWeights = allRoiWeights(TargetRois);
        
        TargetCells = unique([HoloTargets{:}]);
        TargetCells(isnan(TargetCells))=[];
        cellList = 1:numCells;
        notTargetCells = ~ismember(cellList,TargetCells);
        
        TargetCellsWithNan = [HoloTargets{:}];
        
        disp([num2str(numel(TargetCells)) ' Detected Targeted Cells out of ' num2str(numel(TargetRois)) ' ROIs. ' num2str(numel(TargetCells)/numel(TargetRois)*100) '%']);
    catch
        disp('problem with rois, maybe no holo epoch? line 232')
    end
    
    %%stim Params
    numConstant = 1; %Several stim types are shot many times, whereas others are only shot a few times. How many are shot many times
    
    [a b] = sort(uStimCount);
    constIDs = uniqueStims(b(end-numConstant+1:end));
    
    try
        bigList =holoRequests.bigListofSequences;
    catch
        disp('No bigListofSequences, initializing to no rois')
        bigList=[];
    end
    
    uniqueROIs=[];
    for i=1:numel(bigList);
        try
            test = bigList{i}; %changed from Exp Struct to holoRequests.biglist b/c a future epoch overwrites this
        catch
            disp('No bigListofSequences, initializing to no rois')
            test=[];
        end
        uniqueROIs{i} = unique(test);
        
    end
    
    stimParam=[];
    stimNames =[];
    for i =1:numel(uniqueStims)
        ID = uniqueStims(i);
        if ID==0 || ID >=1000
        else
            %         stimNames{i} = ExpStruct.output_names{ID};
            %                 stimParam.Seq(i) = sum(diff(ExpStruct.output_patterns{ID}(:,9))>0);
            %         stimParam.numPulse(i) =  sum(diff(ExpStruct.output_patterns{ID}(:,5))>0);
            stimParam.Seq(i) = sum(diff(ExpStruct.stimlog{ID}{1}(:,7))>0);%formerly 9
            stimParam.numPulse(i) =  sum(diff(ExpStruct.stimlog{ID}{1}(:,5))>0);
            if stimParam.Seq(i)~=0
                try 
                    stimParam.roi{i} = uniqueROIs{stimParam.Seq(i)};
                catch
                    disp('*********Error in stim Params***********')
                    stimParam.roi{i} = nan; 
                end
            else
                stimParam.roi{i} = 0;
            end
        end
    end
    
    %%Determine Vis Stim Trials
    disp('Determine Vis stuff')
    swp =ExpStruct.EpochEnterSweep{DAQepoch};
    if numel(ExpStruct.EpochEnterSweep)<DAQepoch+1
        swp2 = ExpStruct.sweep_counter-1;
    else
        swp2 = ExpStruct.EpochEnterSweep{DAQepoch+1}-1;
    end
    sws = swp:swp2;
    
    if ~exist('DIGITSWEEP')
        DIGITSWEEP = ExpStruct.digitalSweeps;
    end
    
    temp = cellfun(@(x) numel(find(diff(x(:,4))==1)),DIGITSWEEP,'uniformoutput',1);
    visID = temp(sws);
    
    
    temp = cellfun(@(x) find(diff(x(:,6))==1,1),DIGITSWEEP,'uniformoutput',0);
    temp = temp(sws);
    temp = [temp{:}];
    
    timeOfVisStim = temp/20000; %Frame Rate set to daq standard
    visStart = mean(timeOfVisStim(2:end));
    
    temp = cellfun(@(x) find(diff(x(:,6))==-1,1),DIGITSWEEP,'uniformoutput',0);
    temp = temp(sws);
    temp = [temp{:}];
    
    timeOfVisStimEnd = temp/20000; %Frame Rate set to daq standard
    visStop = nanmean(timeOfVisStimEnd(2:end));
    
    uniqueVisStim = unique(visID);
    
    % disp('I Beep'); beep
    disp('Extracting Outs')
    
    autoDataCollector;
end

disp('Saving...')
save([basePath info.date '_' info.mouse '_outfile'], 'out','-v7.3')
disp('All data saved!')

%%

clear mov1 mov2 mov3 sweeps DIGITSWEEP
savePath = basePath;
saveName =['Analysis_e' num2str(DAQepoch) '.mat'];
save(fullfile(savePath,saveName),'-regexp', '^(?!(ExpStruct)$).','-v7.3');
disp('Saved')
