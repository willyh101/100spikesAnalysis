%% Load Experiments/Setup
clear
close all


addpath(genpath('100spikesAnalysis'), genpath('Ian Code'), genpath('analysis-code/matlab')) % will pathing
%%

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

%%
% loadList = loadList(15);


% allLoadList;
% oriLoadList;
% SSTOriLoadList;
% PVOriLoadList;
% u19LoadList;
manifoldLoadList;

% loadPath = 'U:\ioldenburg\outputdata1'
% loadPath = 'C:\Users\ian\Dropbox\Adesnik\Data\outputdata1'
% loadPath = 'C:\Users\SabatiniLab\Dropbox\Adesnik\Data\outputdata1' %Ian Desktop
% loadPath = 'C:\Users\Will\Local Data\100spikes-results\outfiles-ori'
% loadPath = 'T:\Outfiles';
loadPath = 'E:\100spikes-results\outfiles-master';
%%
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

%% Scary hack to overwrite zdf with df data, use with caution!
% you must run allLoadListErrorFixer first because this will mess up
% indexing into the loadList
%%%% DO NOT RUN %%%%%
in =input('Are you sure you want to overwrite (y/n)','s');
if strcmp(in,'y')
All = intentionallyOverwriteZDFWithDF(All);
numExps = numel(All);
else
    disp('Did not run')
end
%%%%%%%%%%%%%%%%%%%%%

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
opts.IndsVisThreshold =-0.20; %default 0.05

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
    & highVisPercentInd ...
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
    ... & numCellsEachEns==10 ...
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

%% Optional group ensembles into small medium and large
numCellsEachEns = outVars.numCellsEachEnsBackup;
numCellsEachEns(numCellsEachEns < 10) = 3;
numCellsEachEns(numCellsEachEns > 10) = 20;

outVars.numCellsEachEns= numCellsEachEns;

%% set them all the same
numCellsEachEns(numCellsEachEns >0 ) = 1;
outVars.numCellsEachEns= numCellsEachEns;

%% Optional Reset ensemble grouping to default
numCellsEachEns = outVars.numCellsEachEnsBackup;
outVars.numCellsEachEns= numCellsEachEns;

%% Basic Response Plots
outVars.defaultColorMap = 'viridis';
plotAllEnsResponse(outVars)
plotResponseBySize(outVars,0)
plotPopResponseBySession(All,outVars)
plotPopResponseByExpressionType(All,outVars);
[All, outVars] = createTSPlotByEnsSize(All,outVars);
[All, outVars] = createTSPlotByEnsSizeAllVis(All,outVars);

%% Distance Response Plots
opts.distBins = 0:25:1000;
plotResponseByDistance(outVars,opts);


%%
ylabel('Pop Response (Mean \Delta Z-dF/F)')
xlabel('Distance From Nearest Target (\mum)')
legend off
%% Second more flexible way to make Distance Plots


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

%% Distance by Vis Response (will only work if consistent number of unique(visID)
plotResponseByDistanceContrast(outVars,opts); %warning won't throw an error even if you have no contrasts
%% Contrast Response Functions

opts.distBinCRF = [50:50:350]; %split Contrast Response Fun by distance
opts.visAlphaCRF = 10.05; %visAlpha for looking just at vis responsive cells;

[outVars] = plotContrastResponseFunction(All,outVars,opts);

%% Orientation Tuning and OSI

[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
[All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically
%% plot vis things
opts.ensOSImethod = 'ensOSI';% 'ensOSI'; 'meanEnsOSI'

plotOSIdists(outVars, opts);
plotPopResponseEnsOSI(outVars, opts)

%% within pyramids cell analysis

opts.visAlpha = 0.05;
[outVars] = makeMeanRespEnsByCell(All, outVars);
opts.restrictToHighOSICells =0.5; 0.5; %0 or 0.5 typical; threshold to restrict tuning analysis to high OSI cells the number is the OSI threshold. set 0 or negative for no restrict
[All, outVars] = compareAllCellsTuning(All, outVars, opts);
[outVars] = getRespByTuningDiff(All, outVars, opts);



%% plot pyr cells connectivity

opts.ensXaxis = 'osi'; % order, osi, dist, corr, size...
plotCompareAllCellsTuning(outVars, opts);
opts.goodOSIthresh = 0.5; %ensemble OSI threshold
plotResponseByDifferenceinAnglePref(outVars, opts)

%% ensemble Plots to look at


%% seperate by posNeg
[All outVars] = posNegIdentifiers(All,outVars,opts);
outVars = getRespByTuningDiffPosNeg(All, outVars);
%%
plotResponseByDiffAnglePrefPosNeg(outVars,opts);

%% ensemble tuning curve examples

OSImin = 0.5;

goodOSIensIdxs = find(outVars.ensOSI > OSImin);


goodOSIens = outVars.ensCurve(2:9,:);
goodOSIensSEM = outVars.ensCurveSEM(2:9,:);
r=5;
% r=505;
% r=573;
% r = randi(numel(goodOSIensIdxs));
% r = goodOSIensIdxs(r);
figure(222)
e = errorbar(goodOSIens(:,r), goodOSIensSEM(:,r));
e.LineWidth = 3;
% e.Color = 'blue';
ylabel('ZDF')
xlabel('Orientation')
xticklabels([0:45:315])
xlim([0.5 8.5])
box off
title('Ensemble Tuning Curve')


ind = outVars.ensIndNumber(r);
hNum = outVars.ensHNumber(r);
% ind=1;
% hNum = 2;
h = All(ind).out.exp.stimParams.roi{hNum};
tg = All(ind).out.exp.holoTargets{h};
tg(isnan(tg))=[];

figure(1222)
clf
nTg = numel(tg);
% [w, h] = splitPretty(nTg, 5, 5);
w=1;h=3;
for i=1:numel(tg)
    c = tg(i);
    subplot(w, h, i)
    curve = All(ind).out.anal.oriCurve(2:9,c);
    curveSEM = All(ind).out.anal.oriCurveSEM(2:9,c);
    e = errorbar(curve, curveSEM);
    e.LineWidth = 2;
    box off
    xlim([0.5 8.5])
    xticklabels([])
    yticklabels([])
    axis off
    title(['Cell ' num2str(i)])
end

%% Red Cell Analysis (will only run if you have the red section on all your recordings).
opts.numExamples = 3;
opts.osiThreshold4Examples = 0.5;
opts.visAlpha = 0.05;
opts.redCellName = 'SST Cells';
[outVars] = plotResponseOfRedCells(All, outVars, opts);
[All, outVars] = redCellTuningAnalysis(All, outVars, opts);
[outVars] = makeMeanRespEns2(All, outVars);

for i=1:numel(outVars.ensStimScore)
    ind = outVars.ensIndNumber(i);
    outVars.isRedByEns{i} = All(ind).out.red.isRed;
end

%% Within Red Cell Analysis

[All, outVars] = compareRedCellsVisResp(All, outVars);
[All, outVars] = compareRedCellsTuning(All, outVars);

opts.redCellName = 'SST';
opts.redCellXaxis = 'size'; % order, osi, dist, corr, size... % correlation stuff is below, so might need to run that first if choosing corr
plotCompareRedCellVisResp(outVars, opts);
plotCompareRedCellVisTuning(outVars, opts);

%% Red cell connectivity
[outVars] = getRespByTuningDiffNotRed(All, outVars);
[outVars] = getRespByTuningDiffRed(All, outVars);
%%
opts.goodOSIthresh = 0.5;
plotResponseByDifferenceinAnglePrefRed(outVars, opts)
%%
opts.goodOSIthresh = 0;
opts.redCellName = 'SST Cells';
opts.goodOSIthresh = 0.001;

plotRedandPyrConnTogether(outVars, opts)

%% Red Distance section
for ind = 1:numExps
    if isfield(All(ind).out.red,'isRed')
        All(ind).out.red.notRed = ~All(ind).out.red.isRed;
    end
end

opts.distType = 'min';

CellToUseVar = 'red.isRed';
[popRespDistEnsRed] = popDistMaker(opts,All,CellToUseVar,0);
CellToUseVar = 'red.notRed';
[popRespDistEnsNotRed] = popDistMaker(opts,All,CellToUseVar,0);

figure(9);clf
ax =subplot(1,2,1);
plotDistRespGeneric(popRespDistEnsRed,outVars,opts,ax);
title('Red Cells')
ax2 = subplot(1,2,2);
plotDistRespGeneric(popRespDistEnsNotRed,outVars,opts,ax2);
title('Not Red Cells')
linkaxes([ax ax2])

%%
figure(42);clf
subplot(1,2,1)
scatter(outVars.ensMaxD,outVars.allRedEnsResp)
refline(0)
subplot(1,2,2)
scatter(outVars.ensMaxD,outVars.popResponseEns)
refline(0)

answerDanQuestionX = outVars.ensMaxD;
answerDanQuestionY = outVars.allRedEnsResp;

% fit(outVars.ensMaxD',outVars.allRedEnsResp','poly1')
%% Pos vs Neg by Distance
opts.posNegThreshold = 0.1;
[All outVars] = posNegIdentifiers(All,outVars,opts);
opts.distType = 'min';
opts.distBins = [0:25:1000];

numEns = numel(outVars.ensStimScore);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlotPos popToPlotNeg
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
    if outVars.ensemblesToUse(i)
    cellToUseVar = outVars.posCellbyInd{i};
    popToPlotPos(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

     cellToUseVar = outVars.negCellbyInd{i};
    popToPlotNeg(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlotPos(i,:) = nan([numel(opts.distBins)-1 1]);
        popToPlotNeg(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')

%%Plot Dist REsp
figure(10);clf
opts.distAxisRang = [0 350];

ax =subplot(1,2,1);
plotDistRespGeneric(popToPlotPos,outVars,opts,ax);
title('Cells that go up')
ax2 =subplot(1,2,2);
plotDistRespGeneric(popToPlotNeg,outVars,opts,ax2);
title('Cells That go down')
% linkaxes([ax ax2]);
% xlim([0 250])
%% Distance of Ensemble
[All, outVars] = defineDistanceTypes(All, outVars);
plotEnsembleDistanceResponse(outVars,100,1)
%% Correlation Pick One. Option A. Vis activity from interleaved Trials
%Functions are Mutually Exclusive.
[All, outVars] = defineCorrelationTypes(All,outVars); %Caution this and below are mutually exclusive
plotEnsembleCorrelationResponse(outVars,200,1);

opts.CorrSpace = linspace(-0.5,0.5,40);
opts.CorrToPlot = 'SpontCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'
[outVars] = plotCorrelationResponse(All,outVars,opts);
%% Correlation Pick One. Option B. Vis Activity from out.vis epoch
[All, outVars] = defineCorrelationTypesOnVis(All, outVars); %Caution this and above are mutually exclusive
plotEnsembleCorrelationResponse(outVars,300,1)

opts.CorrSpace = linspace(-0.2,0.2,40);
opts.CorrToPlot = 'AllCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'
[outVars] = plotCorrelationResponse(All,outVars,opts);

%% Seperate Correlation by Pos and Neg

opts.CorrSpace = linspace(-0.5,0.5,40);linspace(-1,1,40);%
opts.CorrToPlot = 'NoiseCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'

numEns = numel(outVars.posCellbyInd);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlotPos popToPlotNeg
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
    if outVars.ensemblesToUse(i)
    cellToUseVar = outVars.posCellbyInd{i} ;
    popToPlotCorrPos(i,:) = popCorrMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

     cellToUseVar = outVars.negCellbyInd{i};
    popToPlotCorrNeg(i,:) = popCorrMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlotCorrPos(i,:) = nan([numel(opts.CorrSpace)-1 1]);
        popToPlotCorrNeg(i,:) = nan([numel(opts.CorrSpace)-1 1]);
    end
end
disp('Done')

figure(11);clf
ax =subplot(1,2,1);
plotRespGeneric(popToPlotCorrPos,opts.CorrSpace,outVars,opts,ax);
title('Cells that go up')
ax2 =subplot(1,2,2);
plotRespGeneric(popToPlotCorrNeg,opts.CorrSpace,outVars,opts,ax2);
title('Cells That go down')

%% 2D Ensemble Response Plot
opts.distType = 'min';
[outVars] = grandDistanceMaker(opts,All,outVars);
numEns = numel(outVars.ensStimScore);

distRange = [0 35];[50 150];%[50 150];
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

%%

opts.posNegThreshold = 0;0.1;
[All outVars] = posNegIdentifiers(All,outVars,opts);
%% Plot Distance Curves in different Ori Bands
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
% p2 = signrank([datToPlot(:,1);datToPlot(:,5)],datToPlot(:,3));
% [~,p2] = ttest2(datToPlot(:,1),datToPlot(:,3));
% [~,p2] = ttest2([datToPlot(:,1);datToPlot(:,5)],datToPlot(:,3));
title({['Pvalue ' num2str(p)]...
    ['Iso vs Ortho PValue: ' num2str(p2)] })

% title('Cells by Tuning')
% ax2 =subplot(1,2,2);
% plotDistRespGeneric(popToPlotNeg,outVars,opts,ax2);
% title('Cells That go down')


%% ---CLOSE VS FAR---%%
figure(61)
clf

closeDat = squeeze(popToPlot(ensemblesToUse,1,:));

% farbins = find(opts.distBins > 49 & opts.distBins < 101);
farbins = find(opts.distBins > 100);
farbinstart = farbins(1);
farbinstop = farbins(end)-1;
farDat = squeeze(nanmean(popToPlot(ensemblesToUse,farbinstart:farbinstop,:),2));

% err = bootci(10000, {@nanmean, closeDat});
% errorbar(1:numel(m),m,err(1,:),err(2,:))

err = nanstd(closeDat,[],1)./sqrt(size(closeDat,1));
m = nanmean(closeDat);
e1 = errorbar(m,err);
e1.LineWidth = 2;
xlabel('\Delta Preferred Orientation')
ylabel('z-score \DeltaF/F')


hold on

err = nanstd(farDat,[],1)./sqrt(size(farDat,1));
m = nanmean(farDat);
e2 = errorbar(m,err);
e2.LineWidth = 2;


% allDat = squeeze(nanmean(popToPlot(ensemblesToUse,:,:),2));
% err = nanstd(allDat,[],1)./sqrt(size(allDat,1));
% m = nanmean(allDat);
% e3 = errorbar(m,err);
% e3.LineWidth = 2;

% legend({'Close','Far', 'All'})
legend({'Close','Far'})
xticks(1:5)
xticklabels(0:45:180)
box off
ax = gca();
ax.FontSize = 14;

for i=1:5
    [~, p(i)] = ttest2(closeDat(:,i),farDat(:,i));
end
p = p.*5;

    

%% ---ISO VS ORTHO PLOT---%%
figure(62)
clf

% isoClose = mean([closeDat(:,1); closeDat(:,5)],2,'omitnan');
% orthClose = mean(closeDat(:,3),'omitnan');
% 
% isoFar = mean([farDat(:,1); farDat(:,5)],'omitnan');
% orthFar = mean(farDat(:,3),'omitnan');

% isoClose = [closeDat(:,1); closeDat(:,5)];
isoClose = closeDat(:,1);
orthClose = closeDat(:,3);

% isoFar = [farDat(:,1); farDat(:,5)];
isoFar = farDat(:,1);
orthFar = farDat(:,3);

subplot(1,3,1)
clrs = {rgb('RoyalBlue'), rgb('tomato')};
p = plotSpread({isoClose,orthClose},'xnames', {'iso', 'ortho'}, 'showMM',4, 'distributionColors', clrs);
set(findall(gcf(),'type','line'),'markerSize',16)
% set(findall(gcf(),'type','line'),'markerFaceAlpha',0.5)
p{2}(1).Color = 'k';
p{2}(2).Color = 'k';
p{2}(1).LineWidth = 2;
p{2}(2).LineWidth = 2;
uistack(p{2},'top')
title('Close')

[~,p] = ttest2(isoClose,orthClose);
disp(['ttest for iso ortho for close bins: ' num2str(p)])

subplot(1,3,2)
clrs = [rgb('RoyalBlue'); rgb('tomato')];
p = plotSpread({isoFar,orthFar},'xnames', {'iso', 'ortho'}, 'showMM',4, 'distributionColors', clrs);
set(findall(gcf(),'type','line'),'markerSize',16)
% set(findall(gcf(),'type','line'),'markerFaceAlpha',0.5)
p{2}(1).Color = 'k';
p{2}(2).Color = 'k';
p{2}(1).LineWidth = 2;
p{2}(2).LineWidth = 2;
uistack(p{2},'top')
title('Far')

[~,p] = ttest2(isoFar,orthFar);
disp(['ttest for iso ortho for far bins: ' num2str(p)])

subplot(1,3,3)
clrs = [rgb('RoyalBlue'); rgb('tomato')];
% isoAll = reshape(squeeze(nanmean(popToPlot(ensemblesToUse,:,[1,5]),2)),[],1);
isoAll = squeeze(nanmean(popToPlot(ensemblesToUse,:,1),2));
orthAll = squeeze(nanmean(popToPlot(ensemblesToUse,:,3),2));
p = plotSpread({isoAll,orthAll},'xnames', {'iso', 'ortho'}, 'showMM',4, 'distributionColors', clrs);
set(findall(gcf(),'type','line'),'markerSize',16)
% set(findall(gcf(),'type','line'),'markerFaceAlpha',0.5)
p{2}(1).Color = 'k';
p{2}(2).Color = 'k';
p{2}(1).LineWidth = 2;
p{2}(2).LineWidth = 2;
uistack(p{2},'top')
title('All')

[~,p] = ttest2(isoAll,orthAll);
disp(['ttest for iso ortho for all bins: ' num2str(p)])
disp(' ')


%% Plot Distance Plots by criteria
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
%% Two Criteria Distance Plots
%% Plot Distance Plots by criteria
opts.distType = 'min';
opts.distBins =[15:20:150];% [0:25:400];

%things to hold constant
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 ;%& outVars.meanEnsOSI>0.25;
criteria =  outVars.ensMaxD; outVars.ensMeaD;%   outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse);
bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins = [0 400 500 inf];
% bins = [0 400  inf];
% bins = [0 400 500 inf];bins = [0 200 250 inf];

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
% figure(11);
% xlim([0 opts.distBins(end)])
% ylim([-0.1 0.3])
%% OSI stuff
% click to see more examples in subplot 2
ind = 1;
stimNum = 3; % for that expt

% get no stim OSI first
uVisID = unique(All(ind).out.exp.visID);
uVisID(uVisID==0)=[];

uStimID = unique(All(ind).out.exp.stimID);
hNum = uStimID(1);

oriCurveNoStim=[];
oriCurveNoStimSEM = [];
oriTrialNoStimCount =[];
osiNoStim = [];


hNum = uStimID(stimNum);
ensID = find(outVars.ensHNumber == stimNum & outVars.ensIndNumber == ind);
ensPO = outVars.ensPO(ensID);
ensOSI = outVars.ensOSI(ensID);

h = All(ind).out.exp.stimParams.roi{stimNum};
tg = All(ind).out.exp.holoTargets{h};
tg(isnan(tg))=[];
cellList = 1:numel(All(ind).out.anal.ROIinArtifact);

cells2use = All(ind).out.anal.pVisR < 0.05 &...
            ~All(ind).out.anal.ROIinArtifact' &...
            ~All(ind).out.anal.offTargetRisk(h,:) &...
            ~ismember(cellList,tg);

controlStim = uStimID(1);
for i=1:numel(uVisID)
    v= uVisID(i);

    trialsToUse = All(ind).out.exp.visID==v &...
        All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimID==controlStim;

%     disp(sum(trialsToUse))
%
    oriCurveNoStim(i,:) = mean(All(ind).out.exp.rdData(cells2use,trialsToUse), 2);
    oriCurveNoStimSEM(i,:) = sem2(All(ind).out.exp.rdData(cells2use,trialsToUse), 2);
    oriTrialNoStimCount(i) = sum(trialsToUse);
end

oriCurveNoStim = oriCurveNoStim(2:end,:);

oriCurveNoStim = oriCurveNoStim - min(oriCurveNoStim);

% get PO
[~, prefOriNoStim]= max(oriCurveNoStim);

% get OOs
orthoOriNoStim = prefOriNoStim-2;
orthoOriNoStim(orthoOriNoStim<1)=orthoOriNoStim(orthoOriNoStim<1)+8;
orthoOri2NoStim = orthoOriNoStim+4;
orthoOri2NoStim(orthoOri2NoStim>8) = orthoOri2NoStim(orthoOri2NoStim>8)-8;
orthoOriNoStim = cat(1,orthoOriNoStim, orthoOri2NoStim);

for i=1:numel(prefOriNoStim)
    osiNoStim(i) = (oriCurveNoStim(prefOriNoStim(i),i) - mean(oriCurveNoStim(orthoOriNoStim(:,i),i))) / ...
                (oriCurveNoStim(prefOriNoStim(i),i) + mean(oriCurveNoStim(orthoOriNoStim(:,i),i)));
end

% do with stim
oriCurveStim=[];
oriCurveStimSEM = [];
oriTrialStimCount =[];
osiStim = [];


for i=1:numel(uVisID)
    v= uVisID(i);

    trialsToUse = All(ind).out.exp.visID==v &...
        All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimID==hNum;

%     disp(sum(trialsToUse))
%
    oriCurveStim(i,:) = mean(All(ind).out.exp.rdData(cells2use,trialsToUse), 2);
    oriCurveStimSEM(i,:) = sem2(All(ind).out.exp.rdData(cells2use,trialsToUse), 2);
    oriTrialStimCount(i) = sum(trialsToUse);
end

oriCurveStim = oriCurveStim(2:end,:);

oriCurveStim = oriCurveStim - min(oriCurveStim);

% get PO
[~, prefOriStim]= max(oriCurveStim);

% get OOs
orthoOriStim = prefOriStim-2;
orthoOriStim(orthoOriStim<1)=orthoOriStim(orthoOriStim<1)+8;
orthoOri2Stim = orthoOriStim+4;
orthoOri2Stim(orthoOri2Stim>8) = orthoOri2Stim(orthoOri2Stim>8)-8;
orthoOriStim = cat(1,orthoOriStim, orthoOri2Stim);

for i=1:numel(prefOriNoStim)
    osiStim(i) = (oriCurveStim(prefOriStim(i),i) - mean(oriCurveStim(orthoOriStim(:,i),i))) / ...
                (oriCurveStim(prefOriStim(i),i) + mean(oriCurveStim(orthoOriStim(:,i),i)));
end

ensPO = find(ensPO == 0:45:315);
coTunedCellsToEns = prefOriStim == ensPO;


figure(4141)
clf


% ensemble TC
subplot(1,3,1)
e = errorbar(outVars.ensCurve(2:9,ensID), outVars.ensCurveSEM(2:9,ensID));
e.LineWidth = 1.5;
e.Color = 'k';
ylabel('ZDF')
xlabel('Orientation')
% xticklabels(0:45:315)
title('Stimulated Ensemble')
xlim([1 8])


% rand co-tuned cell TC
subplot(1,3,2)
% cellIDs = find(coTunedCellsToEns);
cellIDs = find(osiNoStim > 0.5);
r = randi(numel(cellIDs));
c = cellIDs(r);
e = errorbar(oriCurveNoStim(:,c), oriCurveNoStim(:,c));
e.LineWidth = 1.5;
e.Color = 'blue';
hold on
e2 = errorbar(oriCurveStim(:,c), oriCurveNoStim(:,c));
e2.LineWidth = 1.5;
e2.Color = 'red';
ylabel('ZDF')
xlabel('Orientation')
% xticklabels(0:45:315)
title('Cell Co-Tuned To Ensemble')
legend('No Stim', 'Stim')
xlim([1 8])



% OSI comparison
subplot(1,3,3)
scatter(osiNoStim, osiStim, 'filled')
hold on
scatter(osiNoStim(coTunedCellsToEns), osiStim(coTunedCellsToEns), 'filled')
ylabel('OSI')
xlabel('OSI with stim')
refline(1,0)

%% 2D Plot Maker Corr vs Distance

ensemblesToUse = outVars.ensemblesToUse;
% ensemblesToUse = outVars.ensemblesToUse & outVars.ensOSI>0.66;
ensemblesToUse = outVars.ensemblesToUse & outVars.ensAlCo>0.03;

opts.CorrSpace =linspace(-0.2,0.3,21); %linspace(-0.2,0.3,41);
opts.CorrToPlot = 'AllCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'

opts.distType = 'min';
opts.distBins = [10:20:150];%linspace(0,600,41);% linspace(0,600,91);% [0:25:600];

numEns = numel(outVars.posCellbyInd);
ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popResp2D numResponders
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
    if ensemblesToUse(i)
        cellToUseVar = [];outVars.negCellbyInd{i} ;%outVars.posCellbyInd{i} ;
        [popResp2D(i,:,:) numResponders(i,:,:)] = distCorr2DSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
    else
        popResp2D(i,:,:) = nan([numel(opts.CorrSpace)-1 numel(opts.distBins)-1]);
        numResponders(i,:,:) = nan([numel(opts.CorrSpace)-1 numel(opts.distBins)-1]);
    end
end
disp('done')

%%


colormaptouse = 'rdbu';
climTouse = [-0.2 0.2];%[-0.25 0];%[-0.2 0.2]; %color to use, set this or below to [] to skip.
cprctile = [5 95]; %or both to 0 to auto
xSpcing = 3;12;
ySpcing = 2;4;

numCellsEachEns = outVars.numCellsEachEns;


figure(13);clf
subplot(1,2,1)
imagesc(squeeze(nanmean(popResp2D,1)));
colormap(colormaptouse)

if ~isempty(climTouse)
    caxis(climTouse)
end
ylabel('Correlation')
xlabel('Distance \mum')
xticks(1:xSpcing:numel(opts.distBins))
xticklabels(opts.distBins(1:xSpcing:end))
yticks(1:ySpcing:numel(opts.CorrSpace))
yticklabels(opts.CorrSpace(1:ySpcing:end))

subplot(1,2,2)
imagesc(log10(squeeze(nansum(numResponders,1))));
colormap(colormaptouse)
% caxis([-0.1 0.1])
ylabel('Correlation')
xlabel('Distance \mum')
xticks(1:xSpcing:numel(opts.distBins))
xticklabels(opts.distBins(1:xSpcing:end))
yticks(1:ySpcing:numel(opts.CorrSpace))
yticklabels(opts.CorrSpace(1:ySpcing:end))
h = colorbar;
set(get(h,'label'),'string','Log Number of Sources per Pixel');

figure(16);clf
numSzs = numel(unique(numCellsEachEns(ensemblesToUse)));
szs = unique(numCellsEachEns(ensemblesToUse));

clear ax
for i=1:numSzs
    ax(i) = subplot(1,numSzs,i);
    datToPlot = squeeze(nanmean(popResp2D(numCellsEachEns==szs(i),:,:),1));
    imagesc(datToPlot);
    colormap(colormaptouse)
    if ~isempty(climTouse)
        caxis(climTouse)
    elseif ~isempty(cprctile)
        caxis([prctile(datToPlot(:),cprctile(1)) prctile(datToPlot(:),cprctile(2))]);
    end
    ylabel('Correlation')
    xlabel('Distance \mum')
    xticks(1:xSpcing:numel(opts.distBins))
    xticklabels(opts.distBins(1:xSpcing:end))
    yticks(1:ySpcing:numel(opts.CorrSpace))
    yticklabels(opts.CorrSpace(1:ySpcing:end))

    title(['Ensemble of Size ' num2str(szs(i))]);
end
    linkaxes(ax)
h = colorbar;
set(get(h,'label'),'string','Avg. \Delta Z-Score dF/F');

%% Sparsity

plotSparsityBySize(All,outVars)


%% Stim Rate vs Pop Response
 plotPopRespByStimRate(outVars)

 plotPopRespByNumSpikes(outVars)
 %%
 opts.ensemblesToUseSpikePlot = outVars.ensemblesToUse & outVars.numCellsEachEns ==10;
 opts.plotMeansOnly = 1;
 plotPopRespByNumSpikes2(outVars,opts)
 
 
 %% Section to determine holo/vis interaction
[All] = autodetectVisCondionsInEpoch(All);

temp = arrayfun(@(x) numel(x.out.anal.visCode),All(outVars.IndsUsed),'uniformoutput',1);
indsWithVis = outVars.IndsUsed(temp>1);
ensemblesWithVis = ismember(outVars.ensIndNumber,indsWithVis) & outVars.ensemblesToUse & outVars.ensOSI>0.3;
ensemblesWithVisList = find(ensemblesWithVis);
colorList = {rgb('black') rgb('firebrick')};
titleString = {'No Vis' '\Delta0\circ' '\Delta45\circ' '\Delta90\circ' '\Delta135\circ'...
    '\Delta180\circ' '\Delta225\circ' '\Delta270\circ' '\Delta315\circ' };
oriList = [NaN 0:45:315];

ro = numel(ensemblesWithVisList);
co = 9;
figure(5);clf

for i =1:numel(ensemblesWithVisList)
    ens = ensemblesWithVisList(i);

    ind = outVars.ensIndNumber(ens);
    hNum = outVars.ensHNumber(ens);

    thisEnsOSI = outVars.ensOSI(ens);
    thisEnsPO = outVars.ensPO(ens);

    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);

    subplot(ro,co,1+(i-1)*co)

    trialsToUse = All(ind).out.exp.lowMotionTrials & All(ind).out.exp.lowRunTrials;
    cellsToUse = ~All(ind).out.anal.ROIinArtifact'...
        & ~outVars.offTargetRiskEns{ens}...
       & outVars.isVisR{ind}...
       & outVars.osi{ind} >0.3 ...
        ...& ~(oriList(outVars.prefOris{ind})==thisEnsPO)... | oriList(outVars.prefOris{ind})==thisEnsPO+180 | oriList(outVars.prefOris{ind})==thisEnsPO-180)...
        ;

    datToPlot = All(ind).out.exp.zdfData(cellsToUse,:,trialsToUse & All(ind).out.exp.stimID ==us(1) & All(ind).out.exp.visID == vs(1));
    datToPlot = mean(datToPlot,3);
    fillPlot(datToPlot,[],'ci',colorList{1},'none',colorList{1},0.5);
    datToPlot = All(ind).out.exp.zdfData(cellsToUse,:,trialsToUse & All(ind).out.exp.stimID ==us(hNum) & All(ind).out.exp.visID == vs(1));
    datToPlot = mean(datToPlot,3);
    fillPlot(datToPlot,[],'ci',colorList{2},'none',colorList{2},0.5);

    ylabel(['Ens: ' num2str(ens)]);

    for k = 2:numel(vs)
        if ~isnan(thisEnsPO)
        visCond = All(ind).out.anal.visCode(k);
        v = vs(k);


        visOri = oriList(visCond);
        dif = abs(mod(visOri,180)-mod(thisEnsPO,180));

        cellsToUse2 = cellsToUse & mod(oriList(outVars.prefOris{ind}),180)==mod(visOri,180);


        subplot(ro,co,dif/45+2+(i-1)*co)
        datToPlot = All(ind).out.exp.zdfData(cellsToUse2,:,trialsToUse & All(ind).out.exp.stimID ==us(1) & All(ind).out.exp.visID == vs(k));
        datToPlot = mean(datToPlot,3);
        fillPlot(datToPlot,[],'ci',colorList{1},'none',colorList{1},0.5);
        datToPlot = All(ind).out.exp.zdfData(cellsToUse2,:,trialsToUse & All(ind).out.exp.stimID ==us(hNum) & All(ind).out.exp.visID == vs(k));
        datToPlot = mean(datToPlot,3);
        fillPlot(datToPlot,[],'ci',colorList{2},'none',colorList{2},0.5);
        end
    end
    if i == 1;
        for k=1:co
            subplot(ro,co,k)
            title(titleString{k})
        end
    end

drawnow
end
