%% Analysis template with no DAQ File
clear;
date = '210201';
mouse = 'w29_3';%'I138_1';%'I136_1';
epochs = '1_2_3_4_5_6';

% addpath(genpath('C:\Users\Will\Lab Code\Ian Code'))
% basePath = ['E:\Contrast Modulated Ensembles\' mouse '\' date '\'];
basePath = ['C:\Users\ian\Documents\DATA\F\' mouse '\' date '\'];

path = fullfile(basePath,epochs);

baseName = [mouse '_' date];%'I118a.2_180504';
loadList = {['F_' baseName '_plane1_proc'] ['F_' baseName '_plane2_proc'] ['F_' baseName '_plane3_proc'] ['F_' baseName '_plane4_proc']};

nDepthsTotal = 3;4;%Normally 3;
% 
% physfile = fullfile(basePath,[date '_C' '.mat']);
% % physfile = fullfile(basePath,[date(3:end) '_A' '.mat']);
% try
%     load(physfile)
% catch
%     physfile = fullfile(basePath,[date(3:end) '_A' '.mat']);
%     load(physfile)
% end

%% Experiment

s2pEpoch = 5 ;
DAQepoch = 5 ;

%%  Load without phys
BaseLinePeriod = 1000;
numColors=2;

[localFiles, frameCount, MD, FR, minNumFrames, numFrames, numCells,...
    numFiles, allData, allSPData, allCoM, allDepth, allNP, allNPC,...
    roiMasks, T, BL, allDataCorrect, zData, zSPData, backupData,...
    DFFdata, backupSPData, CellstToExclude, allDataNoNP,...
     ArtifactSizeLeft,...
    ArtifactSizeRight, ROIinArtifact, excludeCells, framesUsed, Tused, zData2]...
    = loadCaDataWithoutPhysFile(path,loadList,DAQepoch,s2pEpoch,nDepthsTotal,BaseLinePeriod,numColors);

%% Determine Stim ID

%Manually get orderBackup
orderWOcatch = cellfun(@(x) x(1),orderBackup);

disp(['Holo Comp recorded ' num2str(numel(orderWOcatch)) ' stims. while ' num2str(numel(localFiles)) ' Files detected']);
disp(['Discrepency: ' num2str(numel(localFiles)-numel(orderWOcatch))])
disc = numel(localFiles)-numel(orderWOcatch);

%% Extract artifact values
val=[];
for i= 1:numel(localFiles);
    tArt = tic;
    fprintf(['Processing file ' num2str(i) ' of ' num2str(numel(localFiles)) ' '])
     v = ScanImageTiffReader(localFiles{i}).data();
        v = permute(v,[2 1 3]);
            gg = v(:,:,1:2:end);
            gg = gg(:,:,1:nDepthsTotal:end);
            temp = gg(:,[40:95 420:470],5:12);
            val(i) = mean(temp(:));
            maxVal(i) = max(temp(:));
            disp(['Took: ' num2str(toc(tArt)) 's.'])
end

%%
smVal = smooth(double(maxVal),100,'lowess')';
rVal = double(maxVal)./smVal;
figure(21);plot(rVal)

catchIDs = find(rVal<0.5); %empirically set
disp([num2str(numel(catchIDs)) ' Catches detected. want ' num2str(disc)]);

stimID = nan([1 numel(localFiles)]);
stimID(catchIDs) = 0;

c=0;
for i =1:numel(stimID);
    if isnan(stimID(i))
        c=c+1;
        stimID(i)=orderWOcatch(c);
    end
end

uniqueStims = unique(stimID); 

%% run Vector
disp('No Run Info')
sz = size(allData);
runVector = zeros(sz([3 2]));
%% holoStuff
holoRequests = holoRequest; % out.holoRequest;
[dfData, zdfData] =  computeDFFwithMovingBaseline(allData);
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

if ~exist('numFiles')
    numFiles = numel(frameCount);
end
numStimuli = 1

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
%% Align Targets Shot to Cells


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


%% vis Section 

disp('navigate to the path with the vis files')
visPath = uigetdir;

%%
temp = dir(localFiles{1});

timeOfFirstTrial= temp.date;
timeOfFirstTrial
temp2 = load(fullfile(visPath,'w29_3_155_004.mat'));
temp2.result.starttime

size(temp2.result.conds,2);

visID1 = nan([1 size(temp2.result.stimParams,2)]);
% visID2 = nan([1 size(temp2.result.stimParams,2)]);

visR1 = temp2.result.stimParams(1,:);
% visR2 = temp2.result.stimParams(3,:);

for i=1:size(temp2.result.conds,2)
    testVal = temp2.result.conds(1,i);
    if ~isnan(testVal)
        visID1(visR1==testVal)=i;
    end
end
%%
visID = visID1; 
    

%% 
clear dat v roiMasks gg temp3 PSTHs tTest temp2 test visID2
savePath = basePath;
saveName =['Analysis_e' num2str(DAQepoch) '.mat'];
save(fullfile(savePath,saveName));
disp('Saved')
