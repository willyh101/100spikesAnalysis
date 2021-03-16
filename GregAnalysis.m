%% Start Here

clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('GregAnalysis.m')); 
addpath(genpath(folder));
rmpath(folder)

%%

%which dataset do you wannt to use
% allLoadList;
 oriLoadList;
% SSTOriLoadList;
% PVOriLoadList;

% loadPath = 'T:\Outfiles';

loadPath = '/Users/gregoryhandy/Research_Local/outputdata1';

%% Load
numExps = numel(loadList);
disp(['There are ' num2str(numExps) ' Exps in this LoadList'])
if numExps ~= 0
    clear All
    if ~iscell(loadList)
        numExps=1;
        temp = loadList;
        clear loadList;
        loadList{1} = temp;
    end
    for ind = 1:numExps
        pTime =tic;
        fprintf(['Loading Experiment ' num2str(ind) '...']);
        All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
        fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
    end
else
    disp('Did you press this by accident?')
end

%% Fix Known Errors
%CAUTION! ERRORS WILL OCCUR IF YOU RUN MORE THAN ONCE!
[All] = allLoadListErrorFixer(All,loadList);

%% clean Data, and create fields.

opts.FRDefault=6;
opts.recWinRange = [0.5 1.5];% %from vis Start in s [1.25 2.5];


%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull

%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded


[All, outVars] = cleanData(All,opts);

ensStimScore        = outVars.ensStimScore;
hzEachEns           = outVars.hzEachEns;
numCellsEachEns     = outVars.numCellsEachEns;
numSpikesEachStim   = outVars.numSpikesEachStim;
percentLowRunTrials = outVars.percentLowRunTrials;
numSpikesEachEns    = outVars.numSpikesEachEns;
numSpikesEachCell   = outVars.numSpikesEachCell;

outVars.numCellsEachEnsBackup = outVars.numCellsEachEns;

names=[];
for Ind = 1:numel(All)
    names{Ind}=lower(strrep(All(Ind).out.info.mouse, '_', '.'));
end
outVars.names = names;

%% Make all dataPlots into matrixes of mean responses
%%Determine Vis Responsive and Process Correlation

opts.visAlpha = 0.05;

%oftarget risk params
opts.thisPlaneTolerance = 11.25;%7.5;%1FWHM%10; %in um;% pixels
opts.onePlaneTolerance = 22.5;%15;%2FWHM %20;
opts.distBins =  [0:25:1000]; [0:25:1000];

[All, outVars] = meanMatrixVisandCorr(All,opts,outVars); %one of the main analysis functions

visPercent = outVars.visPercent;
outVars.visPercentFromExp = visPercent;
ensIndNumber =outVars.ensIndNumber;


%% Optional: Calc pVisR from Visual Epoch [CAUTION: OVERWRITES PREVIOUS pVisR]
[All, outVars] = CalcPVisRFromVis(All,opts,outVars);
visPercent = outVars.visPercent;
outVars.visPercentFromVis = visPercent;

%% if there is a red section (will run even if not...)
[outVars] = detectShotRedCells(All,outVars);
ensHasRed = outVars.ensHasRed;

try
arrayfun(@(x) sum(~isnan(x.out.red.RedCells)),All)
arrayfun(@(x) mean(~isnan(x.out.red.RedCells)),All)
catch end
%% Identify the Experiment type for comparison or exclusion
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;
ensExpressionType = indExpressionType(outVars.ensIndNumber);
outVars.ensExpressionType = ensExpressionType;

%% Missed Target Exclusion Criteria
%detects if too many targets were not detected in S2p
opts.FractionMissable = 0.33; %what percent of targets are missable before you exclude the ens
[outVars] = missedTargetDetector(All,outVars,opts);

ensMissedTargetF = outVars.ensMissedTargetF; %Fraction of targets per ensemble Missed
ensMissedTarget = outVars.ensMissedTarget; %Ensemble is unuseable
numMatchedTargets = outVars.numMatchedTargets;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ensemble Exclusion Section %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% main Ensembles to Use section

numTrialsPerEns =[];numTrialsPerEnsTotal=[]; numTrialsNoStimEns=[];
for ind=1:numExps
    us=unique(All(ind).out.exp.stimID);

    for i=1:numel(us)
        trialsToUse = All(ind).out.exp.lowMotionTrials &...
            All(ind).out.exp.lowRunTrials &...
            All(ind).out.exp.stimSuccessTrial &...
            All(ind).out.exp.stimID == us(i) & ...
            (All(ind).out.exp.visID == 1 | All(ind).out.exp.visID == 0); %restrict just to no vis stim conditions

        numTrialsPerEns(end+1)=sum(trialsToUse);
        numTrialsPerEnsTotal(end+1) = sum(All(ind).out.exp.stimID == us(i));

        if i==1
            numTrialsNoStimEns(ind) = sum(trialsToUse);
        end
    end


end
numTrialsPerEns(numSpikesEachStim==0)=[];
numTrialsPerEnsTotal(numSpikesEachStim==0)=[];

%ID inds to be excluded
opts.IndsVisThreshold = 0.20; %default 0.05

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<opts.IndsVisThreshold)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));

%exclude certain expression types:
uniqueExpressionTypes = outVars.uniqueExpressionTypes;
% excludedTypes = {'AAV CamK2'};
excludedTypes = {'AAV CamK2' 'Ai203'};
% excludedTypes = {'Ai203'};
% excludedTypes = {'AAV Tre'};
% excludedTypes = {};

exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

% only include times where rate == numpulses aka the stim period is 1s.
ensembleOneSecond = outVars.numSpikesEachEns./outVars.numCellsEachEns == outVars.hzEachEns;

%spot to add additional Exclusions
excludeInds = ismember(ensIndNumber,[]); %Its possible that the visStimIDs got messed up

%Options
opts.numSpikeToUseRange = [75 150];[1 inf];[80 120];%[0 1001];
opts.ensStimScoreThreshold = 0.5; % default 0.5
opts.numTrialsPerEnsThreshold = 5; % changed from 10 by wh 4/23 for testing stuff

lowBaseLineTrialCount = ismember(ensIndNumber,find(numTrialsNoStimEns<opts.numTrialsPerEnsThreshold));


ensemblesToUse = numSpikesEachEns > opts.numSpikeToUseRange(1) ...
    & numSpikesEachEns < opts.numSpikeToUseRange(2) ...
    ...& highVisPercentInd ...
    & lowRunInds ...
    & ensStimScore > opts.ensStimScoreThreshold ... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    & ~excludeInds ...
    & numTrialsPerEns > opts.numTrialsPerEnsThreshold ... ;%10;%&...
    & ~lowBaseLineTrialCount ...
    & ~ensHasRed ...
    & ~excludeExpressionType ...
    & ~ensMissedTarget ...
    & numMatchedTargets >= 3 ...
    ...& ensembleOneSecond ... %cuts off a lot of the earlier
    & numCellsEachEns==10 ...    
    ;
%& numCellsEachEns>10 ;

indsSub = ensIndNumber(ensemblesToUse);
IndsUsed = unique(ensIndNumber(ensemblesToUse));

sum(ensemblesToUse)

outVars.ensemblesToUse      = ensemblesToUse;
outVars.IndsUsed            = IndsUsed;
outVars.indsSub             = indsSub;
outVars.numTrialsPerEns     = numTrialsPerEns;
outVars.highVisPercentInd    = highVisPercentInd;
outVars.lowRunInds           = lowRunInds;

%%Optional: Where are the losses comming from

disp(['Fraction of Ens correct Size: ' num2str(mean(numSpikesEachEns > opts.numSpikeToUseRange(1) & numSpikesEachEns < opts.numSpikeToUseRange(2)))]);
disp(['Fraction of Ens highVis: ' num2str(mean(highVisPercentInd))]);
disp(['Fraction of Ens lowRun: ' num2str(mean(lowRunInds))]);
disp(['Fraction of Ens high stimScore: ' num2str(mean(ensStimScore>opts.ensStimScoreThreshold))]);
disp(['Fraction of Ens high trial count: ' num2str(mean(numTrialsPerEns>opts.numTrialsPerEnsThreshold))]);
disp(['Fraction of Control Ens high trial count: ' num2str(mean(~lowBaseLineTrialCount))]);
disp(['Fraction of Ens No ''red'' cells shot: ' num2str(mean(~ensHasRed))]);
disp(['Fraction of Ens usable Expression Type: ' num2str(mean(~excludeExpressionType))]);
disp(['Fraction of Ens enough targets detected by s2p: ' num2str(mean(~ensMissedTarget))]);

disp(['Total Fraction of Ens Used: ' num2str(mean(ensemblesToUse))]);
disp([num2str(sum(ensemblesToUse)) ' Ensembles Included'])


%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;% & ...
%         All(ind).out.exp.visID == 1 | ...
%         All(ind).out.exp.visID == 0;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end
%% Create time series plot

[All, outVars] = createTSPlot(All,outVars);
% [All, outVars] = createTSPlotAllVis(All,outVars);

%% Combine Ensemble Sizes 
numCellsEachEns = outVars.numCellsEachEnsBackup;
outVars.numCellsEachEns= numCellsEachEns;

SizeOption = 1; %0 true; 1 small medium and large; 2 all the same value

switch SizeOption
    case 1
        % Optional group ensembles into small medium and large
        numCellsEachEns = outVars.numCellsEachEnsBackup;
        numCellsEachEns(numCellsEachEns < 10) = 5;
        numCellsEachEns(numCellsEachEns > 10) = 20;
        
        outVars.numCellsEachEns= numCellsEachEns;
    case 2
        % set them all the same
        numCellsEachEns(numCellsEachEns >0 ) = 1;
        outVars.numCellsEachEns= numCellsEachEns;
end

%% Basic Response Plots
outVars.defaultColorMap = 'viridis';
plotAllEnsResponse(outVars)
plotResponseBySize(outVars,0)
plotPopResponseBySession(All,outVars)
plotPopResponseByExpressionType(All,outVars);
[All, outVars] = createTSPlotByEnsSize(All,outVars);
% [All, outVars] = createTSPlotByEnsSizeAllVis(All,outVars);

%% Distance Response Plots
opts.distBins = 0:25:1000;
plotResponseByDistance(outVars,opts);


%% relabel for printing
ylabel('Pop Response (Mean \Delta Z-dF/F)')
xlabel('Distance From Nearest Target (\mum)')
legend off

%% Compare Distance responses
figure(102);clf

distTypes = {'min' 'geo' 'mean' 'harm'};
for i =1:4
    disp(['working on ' distTypes{i}])
    opts.distType = distTypes{i}; %options: min geo mean harm
    CellToUseVar = [];
    [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
    ax = subplot(2,2,i);
    opts.distAxisRange = [0 550]; %[0 350] is stand
    plotDistRespGeneric(popRespDist,outVars,opts,ax);
    title(distTypes{i})
end
disp('done')

%% Distance of Ensemble
[All, outVars] = defineDistanceTypes(All, outVars);
plotEnsembleDistanceResponse(outVars,100,1)

%% Older Correlation Based analyses
% %% Correlation Pick One. Option A. Vis activity from interleaved Trials
% %Functions are Mutually Exclusive.
% [All, outVars] = defineCorrelationTypes(All,outVars); %Caution this and below are mutually exclusive
% plotEnsembleCorrelationResponse(outVars,200,1);
% 
% opts.CorrSpace = linspace(-0.5,0.5,40);
% opts.CorrToPlot = 'SpontCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'
% [outVars] = plotCorrelationResponse(All,outVars,opts);

% %% Correlation Pick One. Option B. Vis Activity from out.vis epoch
% [All, outVars] = defineCorrelationTypesOnVis(All, outVars); %Caution this and above are mutually exclusive
% plotEnsembleCorrelationResponse(outVars,300,1)
% 
% opts.CorrSpace = linspace(-0.2,0.2,40);
% opts.CorrToPlot = 'AllCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'
% [outVars] = plotCorrelationResponse(All,outVars,opts);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 2D Plot Maker Corr vs Distance
% 
% ensemblesToUse = outVars.ensemblesToUse;
% % ensemblesToUse = outVars.ensemblesToUse & outVars.ensOSI>0.66;
% ensemblesToUse = outVars.ensemblesToUse & outVars.ensAlCo>0.03;
% 
% opts.CorrSpace =linspace(-0.2,0.3,21); %linspace(-0.2,0.3,41);
% opts.CorrToPlot = 'AllCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'
% 
% opts.distType = 'min';
% opts.distBins = [10:20:150];%linspace(0,600,41);% linspace(0,600,91);% [0:25:600];
% 
% numEns = numel(outVars.posCellbyInd);
% ensIndNumber = outVars.ensIndNumber;
% ensHNumber = outVars.ensHNumber;
% 
% clear popResp2D numResponders
% for i=1:numEns %i know its slow, but All is big so don't parfor it
%     if mod(i,round(numEns/10))==1
%         fprintf('.')
%     end
%     if ensemblesToUse(i)
%         cellToUseVar = [];outVars.negCellbyInd{i} ;%outVars.posCellbyInd{i} ;
%         [popResp2D(i,:,:) numResponders(i,:,:)] = distCorr2DSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
%     else
%         popResp2D(i,:,:) = nan([numel(opts.CorrSpace)-1 numel(opts.distBins)-1]);
%         numResponders(i,:,:) = nan([numel(opts.CorrSpace)-1 numel(opts.distBins)-1]);
%     end
% end
% disp('done')
% 
% %% make Plot, seperate because the generate data section is slow. 
% 
% 
% colormaptouse = 'rdbu';
% climTouse = [-0.2 0.2];%[-0.25 0];%[-0.2 0.2]; %color to use, set this or below to [] to skip.
% cprctile = [5 95]; %or both to 0 to auto
% xSpcing = 3;12;
% ySpcing = 2;4;
% 
% numCellsEachEns = outVars.numCellsEachEns;
% 
% 
% figure(13);clf
% subplot(1,2,1)
% imagesc(squeeze(nanmean(popResp2D,1)));
% colormap(colormaptouse)
% 
% if ~isempty(climTouse)
%     caxis(climTouse)
% end
% ylabel('Correlation')
% xlabel('Distance \mum')
% xticks(1:xSpcing:numel(opts.distBins))
% xticklabels(opts.distBins(1:xSpcing:end))
% yticks(1:ySpcing:numel(opts.CorrSpace))
% yticklabels(opts.CorrSpace(1:ySpcing:end))
% 
% subplot(1,2,2)
% imagesc(log10(squeeze(nansum(numResponders,1))));
% colormap(colormaptouse)
% % caxis([-0.1 0.1])
% ylabel('Correlation')
% xlabel('Distance \mum')
% xticks(1:xSpcing:numel(opts.distBins))
% xticklabels(opts.distBins(1:xSpcing:end))
% yticks(1:ySpcing:numel(opts.CorrSpace))
% yticklabels(opts.CorrSpace(1:ySpcing:end))
% h = colorbar;
% set(get(h,'label'),'string','Log Number of Sources per Pixel');
% 
% figure(16);clf
% numSzs = numel(unique(numCellsEachEns(ensemblesToUse)));
% szs = unique(numCellsEachEns(ensemblesToUse));
% 
% clear ax
% for i=1:numSzs
%     ax(i) = subplot(1,numSzs,i);
%     datToPlot = squeeze(nanmean(popResp2D(numCellsEachEns==szs(i),:,:),1));
%     imagesc(datToPlot);
%     colormap(colormaptouse)
%     if ~isempty(climTouse)
%         caxis(climTouse)
%     elseif ~isempty(cprctile)
%         caxis([prctile(datToPlot(:),cprctile(1)) prctile(datToPlot(:),cprctile(2))]);
%     end
%     ylabel('Correlation')
%     xlabel('Distance \mum')
%     xticks(1:xSpcing:numel(opts.distBins))
%     xticklabels(opts.distBins(1:xSpcing:end))
%     yticks(1:ySpcing:numel(opts.CorrSpace))
%     yticklabels(opts.CorrSpace(1:ySpcing:end))
% 
%     title(['Ensemble of Size ' num2str(szs(i))]);
% end
%     linkaxes(ax)
% h = colorbar;
% set(get(h,'label'),'string','Avg. \Delta Z-Score dF/F');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start of the OSI/Orientation Section %%%
%% Orientation Tuning and OSI %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
[All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically
%% plot vis things
opts.ensOSImethod = 'ensOSI';% 'ensOSI'; 'meanEnsOSI'

plotOSIdists(outVars, opts);
plotPopResponseEnsOSI(outVars, opts)

%% within pyramids cell analysis 
%Compare Population response (regardless of distance) by orientation
%similarity

opts.visAlpha = 0.05;
[outVars] = makeMeanRespEnsByCell(All, outVars);
opts.restrictToHighOSICells =0.5; 0.5; %0 or 0.5 typical; threshold to restrict tuning analysis to high OSI cells the number is the OSI threshold. set 0 or negative for no restrict
[All, outVars] = compareAllCellsTuning(All, outVars, opts);
[outVars] = getRespByTuningDiff(All, outVars, opts);
%%plot pyr cells connectivity

opts.ensXaxis = 'osi'; % order, osi, dist, corr, size...
plotCompareAllCellsTuning(outVars, opts);
opts.goodOSIthresh = 0.5; %ensemble OSI threshold
plotResponseByDifferenceinAnglePref(outVars, opts)

%% 2D Ensemble Response Plot
opts.distType = 'min';
[outVars] = grandDistanceMaker(opts,All,outVars);
numEns = numel(outVars.ensStimScore);

distRange = [0 35]%[50 150];%[50 150];
for i = 1:numEns
    responseToPlot(i) = nanmean(outVars.mRespEns{i}(...
        outVars.distToEnsemble{i}>=distRange(1) & ...
        outVars.distToEnsemble{i}<distRange(2) & ...
        ~outVars.offTargetRiskEns{i} ...
        ...& outVars.isRedByEns{i} ...
        ));
end
% responseToPlot = outVars.popResponseEns; %outVars.mRespEns;
yToPlot = outVars.ensOSI;
xToPlot = outVars.ensMaxD;

ensemblesToUse = outVars.ensemblesToUse & outVars.numMatchedTargets>2 ;%& outVars.numMatchedTargets<12;
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;% &  outVars.ensOSI>-0.25 & outVars.numMatchedTargets>=3 & outVars.ensMaxD<400;% outVars.numMatchedTargets>=3 &    ;

figure(142);clf
scatter(xToPlot(ensemblesToUse),yToPlot(ensemblesToUse),100,responseToPlot(ensemblesToUse),'filled');
colormap(puor)
colorbar
% caxis([-0.1 0.1])
caxis([-0.25 0.25])
xlabel({'Span of Ensemble (\mum)' ' (i.e. max distance between targets)'})
ylabel('Ensemble OSI')

%% interp response of all cells split into different condition bands
%% Plot Sorted Cell Response

ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 & outVars.ensMaxD<400 & outVars.ensOSI>-0.5;% & outVars.numMatchedTargets>=3;%    ;
yToPlot =outVars.ensOSI; outVars.ensMaxD; %outVars.ensOSI; outVars.numCellsEachEns;outVars.ensOSI;

numpts =1000;

distRanges = [0 30; 31 60; 61 inf];

% distRange = [0 inf];[0 inf];[0 50];%[50 150];
sortedCellResp = [];

figure(1245)
clf
for i=1:3
    subplot(1,3,i)
distRange = distRanges(i,:);
for k = 1:numEns
    fprintf('.');
    if mod(k,100)==0;
        disp(' ')
    end
    cellsToUse =         outVars.distToEnsemble{k}>=distRange(1) & ...
        outVars.distToEnsemble{k}<distRange(2) & ...
        ~outVars.offTargetRiskEns{k};
    if sum(cellsToUse)<5
        temp = [NaN NaN];
        
       ensemblesToUse(k)=0;
    else
        temp = outVars.mRespEns{k}(cellsToUse);
        temp = sort(temp);
    end
    sortedCellResp(k,:) = interp1(1:numel(temp),temp,linspace(numel(temp),1,numpts));
end
disp('done')

ensemblesToUseList = find(ensemblesToUse);

[s sidx] = sort(yToPlot(ensemblesToUse)) ;
sortedEnsList = ensemblesToUseList(sidx);

datToPlot = sortedCellResp(sortedEnsList,:);
% figure(143);clf
imagesc(datToPlot')
colormap(rdbu)
caxis([-.4 .4])
title(['Distances ' num2str(distRange)])

xticks([1 numel(s)])
xticklabels({'Close' 'Far'})
yticks([])
ylabel('Cell Response Sorted')
xlabel('Span of Ensemble')

box off
end
% colorbar
axis tight


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Distance Curves by the difference in preffered Ori responder to ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.distType = 'min';
opts.distBins =[10:20:150];% [0:25:400];

%this is where you change the criteria of what ensembles are included
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI>0.66;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI<0.23;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.meanEnsOSI>0.5 &  outVars.ensOSI>0.25  & outVars.ensMaxD<475 ;%& outVars.numMatchedTargets>=3;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI>0.68  & outVars.ensMaxD>500 ;%& outVars.numMatchedTargets>=3;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI>0.66  & outVars.ensMaxD<400;
sum(ensemblesToUse)

oriVals = [NaN 0:45:315];
% numEns = numel(outVars.posCellbyInd);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end

    diffsPossible = [0 45 90 135 180];

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

       cellOris = oriVals(outVars.prefOris{ind});
       cellOrisDiff = abs(cellOris-outVars.ensPO(i));
       cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);


       %%This is where you change the criteria of what cells are included
    cellToUseVar = ~outVars.offTargetRiskEns{i}...
        & outVars.pVisR{ind} < 0.05 ...
        & outVars.osi{ind} > 0.25 ...
        ... & outVars.posCellbyInd{i} ... %if you want to only include cells that went up
        ...& outVars.isRedByEns{i} ...  %if you want to excluded red cells (i.e. interneurons) 
        ;

    for k=1:numel(diffsPossible)
        popToPlot(i,:,k) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar & abs(cellOrisDiff)==diffsPossible(k),0,ensHNumber(i));
%         popToPlot(i,:,k) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar ,0,ensHNumber(i));

    end


    else
        popToPlot(i,:,:) = nan([numel(opts.distBins)-1 1 numel(diffsPossible)]);
    end
end
disp('Done')

%%Plot Dist REsp
figure(10);clf
opts.distAxisRang = [0 350];

figure(11);clf;hold on
% ax =subplot(1,1,1);
%  hold on
colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
for k = 1:numel(diffsPossible)
    figure(10);
    ax(k) =subplot(1,numel(diffsPossible),k);
    title(['Cells Pref Angle \Delta' num2str(diffsPossible(k)) '\circ'])
    [eHandle] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,ax(k));
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end
    ylabel('Pop Response (Mean \DeltaF/F)')

    figure(11); tempax = subplot(1,1,1);
    [eHandle] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,tempax);
    eHandle{1}.Color = colorListOri{k};
        ylabel('Pop Response (Mean \DeltaF/F)')

end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])

figure(10);
xlim([0 opts.distBins(end)])


figure(12);clf
datToPlot = squeeze(popToPlot(ensemblesToUse,1,:));
ciToPlot = nanstd(datToPlot,[],1)./sqrt(size(datToPlot,1))*1.93;
errorbar(nanmean(datToPlot),ciToPlot);
p = anova1(datToPlot,[],'off');
p2 = signrank(datToPlot(:,1),datToPlot(:,3));
title({['Pvalue ' num2str(p)]...
    ['Iso vs Ortho PValue: ' num2str(p2)] })


%% Plot Distance Plots by criteria
% More Criteria in this case set up as ensOSI, but its supposed to be
% flexible in case you want to change it later.

opts.distType = 'min';
opts.distBins =[10:20:150];% [0:25:400];

%things to hold constant
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;%  &  outVars.ensOSI>0.75;;
criteria = outVars.ensOSI;outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse);
bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];

% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI>0.75;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensMaxD<300;% outVars.numMatchedTargets>=3 &    ;

% criteria = outVars.ensOSI; %outVars.ensMaxD;
% useableCriteria = criteria(ensemblesToUse);
% bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

    cellToUseVar = ~outVars.offTargetRiskEns{i}...
        & outVars.pVisR{ind} < 0.1 ...
        & outVars.osi{ind} > 0.25 ...
        ;

        popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')

figure(10);clf
figure(11);clf

colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
for k = 1:numel(bins)-1
    figure(10);ax(k) =subplot(1,numel(bins)-1,k);
    title(['Ens Criteria: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ])
    popToPlotTemp = popToPlot;
    popToPlotTemp(criteria<bins(k) | criteria>bins(k+1),:)=NaN;
    [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(k));
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end

    figure(11); tempax = subplot(1,1,1);
    [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
    eHandle{1}.Color = colorListOri{k};
end
linkaxes(ax)
figure(11);
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])
figure(10);
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Two Criteria Distance Plots
%%% Plot Distance Plots by Ensemble OSI and Distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


opts.distType = 'min';
opts.distBins =[15:20:150];% [0:25:400];

%adjust the Ensembles used 
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 ;%& outVars.meanEnsOSI>0.25;
criteria =  outVars.ensMaxD; outVars.ensMeaD;%   outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse);
bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins = [0 400 500 inf]; 
% bins = [0 200 250 inf]; 

% bins = [linspace(min(useableCriteria),max(useableCriteria),5)];
criteria2 = outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria2(ensemblesToUse);
bins2 = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins2 = [0 0.33 0.7 inf]; 

% bins2 = [linspace(min(useableCriteria),max(useableCriteria),5)];

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

        %adjust the responder cells included here
    cellToUseVar = ~outVars.offTargetRiskEns{i}...
         & outVars.pVisR{ind} < 0.05 ...
        ...& outVars.osi{ind} > 0.25 ...
        ... & outVars.posCellbyInd{i} ...
        ;

        popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')

figure(13);clf

colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
c=0;
for i = 1:numel(bins2)-1
for k = 1:numel(bins)-1
    c=c+1;
    
    ax(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
popToPlotTemp = popToPlot;
ensembleExcluder  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;
popToPlotTemp(~ensembleExcluder,:)=NaN;
    [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(c));
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end
        title({...
            ['Ens Dist: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            ['Ens OSI: ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            ['Num Ens: ' num2str(sum(ensembleExcluder & ensemblesToUse)) ] ...
            } )

    ylabel('Pop Response (Mean \DeltaF/F)')

%     figure(11); tempax = subplot(1,1,1);
%     [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
%     eHandle{1}.Color = colorListOri{k};
end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.2 0.25])