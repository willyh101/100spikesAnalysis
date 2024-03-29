%% Load Experiments/Setup
clear
close all


addpath(genpath('100spikesAnalysis'))
%% loadLists

oriLoadList;
allLoadList;

% loadPath = 'path/to/outfiles/directory';
loadPath = 'T:\Outfiles';

addpath(genpath(loadPath))

%% Load data

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

%% error fixer
%CAUTION! ERRORS WILL OCCUR IF YOU RUN MORE THAN ONCE!
[All] = allLoadListErrorFixer(All,loadList);

%% Set Data To use
for ind=1:numExps
    All(ind).out.exp.dataToUse = All(ind).out.exp.dfData;
end
disp('Data To Use is set')

%% clean Data, and create fields.

opts.FRDefault=6;
opts.recWinRange = [0.5 1.5]; %[0.5 1.5];[1.5 2.5];%[0.5 1.5];% %from vis Start in s [1.25 2.5];


%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull
opts.stimSuccessByZ = 1; %if 0 calculates based on dataTouse, if 1 calculates based on zdfDat;

%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded
opts.runValPercent = 0.75; %percent of frames that need to be below run threshold

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

%% restrict Cells to use
opts.minMeanThreshold = 0.25;
opts.maxMeanThreshold = inf;

opts.verbose =0;
[All, cellExcludeResults] = cellExcluder(All,opts); 
allResults = cat(1,cellExcludeResults{:});
disp(['In total ' num2str(sum(allResults)) ' Cells Excluded. ' num2str(mean(allResults)*100,2) '%']);
disp(['Overall ' num2str(sum(~allResults)) ' Cells Passed!'])

opts.minNumCellsInd=250;
tooFewCellsInds = cellfun(@(x) sum(~x)<opts.minNumCellsInd,cellExcludeResults);
disp([ num2str(sum(tooFewCellsInds)) ' inds have < ' num2str(opts.minNumCellsInd) ' cells, and should be exccluded']);



%% Make all dataPlots into matrixes of mean responses
%%Determine Vis Responsive and Process Correlation

opts.visAlpha = 0.05;

%oftarget risk params
opts.thisPlaneTolerance = 11.25;%7.5;%1FWHM%10; %in um;% pixels
opts.onePlaneTolerance = 22.5;%15;%2FWHM %20;
opts.distBins =  [0:25:1000]; [0:25:1000];
opts.skipVis =1;

[All, outVars] = meanMatrixVisandCorr(All,opts,outVars); %one of the main analysis functions

visPercent = outVars.visPercent;
outVars.visPercentFromExp = visPercent;
ensIndNumber =outVars.ensIndNumber;



%% Misc Additional Variables:
%%RedSection: if there is a red section (will run even if not...)
[outVars] = detectShotRedCells(All,outVars);
ensHasRed = outVars.ensHasRed;

try
arrayfun(@(x) sum(~isnan(x.out.red.RedCells)),All)
arrayfun(@(x) mean(~isnan(x.out.red.RedCells)),All)
catch end
%%Identify the Experiment type for comparison or exclusion
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;
ensExpressionType = indExpressionType(outVars.ensIndNumber);
outVars.ensExpressionType = ensExpressionType;

%%Missed Target Exclusion Criteria
%detects if too many targets were not detected in S2p
opts.FractionMissable = 0.33; %what percent of targets are missable before you exclude the ens
[outVars] = missedTargetDetector(All,outVars,opts);

ensMissedTargetF = outVars.ensMissedTargetF; %Fraction of targets per ensemble Missed
ensMissedTarget = outVars.ensMissedTarget; %Ensemble is unuseable
numMatchedTargets = outVars.numMatchedTargets;
%%Determine Date
ensDate=[];
for i = 1:numel(ensIndNumber)
    ensDate(i) = str2num(All(ensIndNumber(i)).out.info.date);
end
outVars.ensDate=ensDate;

%%Identify duplicate holograms
[outVars] = identifyDuplicateHolos(All,outVars);

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
opts.IndsVisThreshold = 0.05; %default 0.05

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<opts.IndsVisThreshold)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));
lowCellCount = ismember(ensIndNumber,find(tooFewCellsInds));

%exclude certain expression types:
uniqueExpressionTypes = outVars.uniqueExpressionTypes;
excludedTypes ={'AAV CamK2' 'Ai203' 'neo-IV Tre 2s' 'IUE CAG' 'SepW1 CAG 2s' 'control' };


exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

% only include times where rate == numpulses aka the stim period is 1s.
ensembleOneSecond = outVars.numSpikesEachEns./outVars.numCellsEachEns == outVars.hzEachEns;

%spot to add additional Exclusions
excludeInds = ismember(ensIndNumber,[]); %Its possible that the visStimIDs got messed up
% excludeInds = ismember(ensIndNumber,[]); 

%Options
opts.numSpikeToUseRange = [90 110];[1 inf];[80 120];%[0 1001];
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
    & ensembleOneSecond ... %cuts off a lot of the earlier
    ... & (numCellsEachEns==10 | numCellsEachEns==5 | numCellsEachEns==20) ...
    & (numCellsEachEns==10 | numCellsEachEns==3 | numCellsEachEns==33) ...
    ... & (numCellsEachEns==10 | numCellsEachEns==3 | numCellsEachEns==33 | numCellsEachEns==5 | numCellsEachEns==20) ...
    ... & numCellsEachEns==10 ...
    & ensDate >= -210428 ...
    ...& outVars.hzEachEns == 10 ...
    ...& (outVars.hzEachEns >= 9 & outVars.hzEachEns <= 12 ...
    & ~lowCellCount ...
    ;

%%remove repeats
 [ensemblesToUse, outVars] = removeRepeatsFromEnsemblesToUse(ensemblesToUse,outVars);

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
disp(['Fraction of Ens number targets matched >=3: ' num2str(mean(numMatchedTargets >= 3))]);
disp(['Fraction of Ens Stim took 1s (aka correct stim Rate): ' num2str(mean(ensembleOneSecond))]);
disp(['Fraction of Ens that were not repeats: ' num2str(mean(~outVars.removedRepeats)) ]);
disp(['Fraction of Ens high Cell Count: ' num2str(mean(~lowCellCount))]);


disp(['Total Fraction of Ens Used: ' num2str(mean(ensemblesToUse))]);
% disp([num2str(sum(ensemblesToUse)) ' Ensembles Included'])
disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);


%% Ensure Ensembles come from inds where each type of holos was presented.

ensSizes =[3 5 10 20 33]; [3 10 33];[5 10 20];
indsBySize =[];
for i = 1:numel(ensSizes);
    indsBySiz{i} = unique(ensIndNumber(ensemblesToUse & numCellsEachEns ==ensSizes(i)));
end
AllList = 1:numel(All);
temp = cellfun(@(x) ismember(AllList,x),indsBySiz,'uniformoutput',0);
indsWithAllSize = find(all(cat(1,temp{:})));

indsWithAllSize = find(sum(cat(1,temp{:}))>=2);

ensemblesToUse = ensemblesToUse & ismember(ensIndNumber,indsWithAllSize);
disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);

outVars.ensemblesToUse = ensemblesToUse;
outVars.IndsUsed = unique(ensIndNumber(ensemblesToUse));

%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials ...
        & All(ind).out.exp.lowRunTrials ...
        & All(ind).out.exp.stimSuccessTrial ...
        & (All(ind).out.exp.visID == 1 |  All(ind).out.exp.visID == 0 ) ...
            ;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end

%% Create time series plot

[All, outVars] = createTSPlot(All,outVars);
% [All, outVars] = createTSPlotAllVis(All,outVars);

%% Optional group ensembles into small medium and large
numCellsEachEns = outVars.numCellsEachEnsBackup;
numCellsEachEns(numCellsEachEns < 10) = 5;
numCellsEachEns(numCellsEachEns > 12) = 20;

outVars.numCellsEachEns= numCellsEachEns;

%% set them all the same
numCellsEachEns(numCellsEachEns >0 ) = 1;
outVars.numCellsEachEns= numCellsEachEns;

%% Optional Reset ensemble grouping to default
numCellsEachEns = outVars.numCellsEachEnsBackup;
outVars.numCellsEachEns= numCellsEachEns;
%%
muPerPx = 800/512;
opts.thisPlaneTolerance =15/muPerPx;% 15/muPerPx;
opts.onePlaneTolerance = 30/muPerPx; %30/muPerPx;

recalcOffTargetRisk;
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

%% Just a few with different binning
figure(104);clf;
dataInPlots =[];

ax = subplot(1,1,1);
opts.distType = 'min';
opts.distBins = 15:15:250; %can be set variably 0:25:1000 is defaultt
opts.distAxisRange = [0 250]; %[0 350] is stand
% CellToUseVar = 'anal.cellsToInclude & All(ind).out.anal.pVisR<0.05';%[];
CellToUseVar = 'anal.cellsToInclude ';%[];

[popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
[eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
dataInPlots{1}=outDat{1};
for i =1: size(eHandle,2)-1
eHandle{i}.CapSize =0;
end
title('min')

%% Compare Distance responses
figure(102);clf

distTypes = {'min' 'geo' 'mean' 'harm' 'median' 'centroid'};
for i =1:6
    disp(['working on ' distTypes{i}])
    opts.distType = distTypes{i}; %options: min geo mean harm
    CellToUseVar = 'anal.cellsToInclude';
    [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
    ax = subplot(2,3,i);
    opts.distAxisRange = [0 350]; %[0 350] is stand
    plotDistRespGeneric(popRespDist,outVars,opts,ax);
    title(distTypes{i})
    drawnow
end
disp('done')

