function stm = ExportStm(inputParams)

path = inputParams.fullPth;
physfile = inputParams.physfile;
loadList = inputParams.loadList;
DAQepoch = inputParams.DAQepoch;
s2pEpoch = inputParams.s2pEpoch;


BaseLinePeriod = 1000;

if ~exist('ExpStruct')
    disp('Reloading ExpStruct')
    load(physfile,'ExpStruct');
end

if ~exist('DIGITSWEEP')
    DIGITSWEEP = ExpStruct.digitalSweeps;
end

numColors=2;
nDepthsTotal = 3;

disp('The Scary Load Part...')
tScary = tic;
[localFiles, frameCount, MD, FR, minNumFrames, numFrames, numCells,...
    numFiles, allData, allSPData, allCoM, allDepth, allNP, allNPC,...
    roiMasks, T, BL, allDataCorrect, zData, zSPData, backupData,...
    DFFdata, backupSPData, ExpStruct, CellstToExclude, allDataNoNP,...
    stimID, uniqueStims, uStimCount, runVector, ArtifactSizeLeft,...
    ArtifactSizeRight, ROIinArtifact, excludeCells, framesUsed, Tused, zData2]...
    = loadCaData(path,loadList,DAQepoch,s2pEpoch,nDepthsTotal,ExpStruct,DIGITSWEEP,BaseLinePeriod, numColors);
swp =ExpStruct.EpochEnterSweep{DAQepoch};
hr = ExpStruct.Holo.Sweeps_holoRequestNumber(swp);

holoRequests = ExpStruct.Holo.holoRequests{hr};

if ~isfield(holoRequests,'DE_list')
    try
       holoRequests = ExpStruct.Holo.holoRequests{hr-1};
    catch
    end
end

disp('Computing df and zdf...')
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

offsets = -mean(cat(1,Toffsets{:}))
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


    TargetRois = unique(cat(2,roisTargets{:}));
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
                stimParam.roi{i} = 0;
                disp('stimParam Issues')
            end
        else
            stimParam.roi{i} = 0;
        end
    end
end

%%Determine Vis Stim Trials

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

disp('done process. now extracting stm')

%% stimTest
stm.DAQepoch = DAQepoch;
stm.zdfData = zdfData;
stm.allData = allData;
stm.stimID =stimID;
% 
% swpStart = ExpStruct.EpochEnterSweep{DAQepoch};
% Hnum = ExpStruct.Holo.Sweeps_holoRequestNumber(swpStart);
% HR =ExpStruct.Holo.holoRequests{Hnum};

stm.holoRequest = holoRequests; %HR;

stm.runVal = runVector;
stm.lowMotionTrials = lowMotionTrials;

stm.holoTargets = HoloTargets;
stm.rois = roisTargets;
stm.allCoM = allCoM;
stm.allDepth = allDepth;
stm.stimCoM = stimCoM;
stm.stimDepth = stimDepth;
stm.targetedCells = targettedCells;
stm.uniqueStims = uniqueStims;
stm.outputsInfo = outputPatternTranslator(ExpStruct,uniqueStims);

tempOutputOrder = stm.outputsInfo.OutputOrder;
tempOutputOrder(tempOutputOrder==0)=[];
% exp.output_names = ExpStruct.output_names;

stm.stimParams = stimParam;
try
stm.stimParams.Hz = holoRequests.holoStimParams.hzList(tempOutputOrder);
stm.stimParams.numCells = holoRequests.holoStimParams.cellsPerHolo(tempOutputOrder);
stm.stimParams.powers = holoRequests.holoStimParams.powerList(tempOutputOrder);
catch
    disp('no holoStimParams')
end

stm.Tarray = Tarray; %Motion Correct trace;
stm.dfData = dfData; %non zscored data;

stm.offsets = offsets; %sometimes different epochs calc subtly different offsets, recorded here
