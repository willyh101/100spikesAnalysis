%% Load Experiments/Setup
clear
close all

masterTic = tic;

addpath(genpath('100spikesAnalysis'))
%% loadLists

% oriLoadList;
sepW1LoadList;
% % allLoadList;

% loadPath = 'path/to/outfiles/directory';
loadPath = 'T:\Outfiles';
loadPath = 'C:\Users\ian\Dropbox\Outfiles'

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


%% Mark and exclude files that don't have opsin negative data
noOpsinDataInd=[];
for ind =1:numExps
   if ~isfield(All(ind).out.info, 'value1020')
        noOpsinDataInd = [noOpsinDataInd ind];
   end
end


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


%% REQUIRED: Calc pVisR from Visual Epoch [CAUTION: OVERWRITES PREVIOUS pVisR]
% Always do this!! not all experiments had full orientation data during the
% experiment epoch (but did during the vis epoch)
disp('Recalculating Vis stuff...')
opts.visRecWinRange = [0.5 1.5]; [0.5 1.5];
[All, outVars] = CalcPVisRFromVis(All,opts,outVars);
visPercent = outVars.visPercent;
outVars.visPercentFromVis = visPercent;

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
excludedTypes ={'AAV CamK2' 'Ai203' 'neo-IV Tre 2s' 'IUE CAG'}; %'SepW1 CAG 2s'


exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

% only include times where rate == numpulses aka the stim period is 1s.
ensembleOneSecond = outVars.numSpikesEachEns./outVars.numCellsEachEns == outVars.hzEachEns;

%spot to add additional Exclusions
excludeInds = ismember(ensIndNumber,noOpsinDataInd); %Its possible that the visStimIDs got messed up
% excludeInds = ismember(ensIndNumber,[]); 

%Options
opts.numSpikeToUseRange = [1 inf];[80 120];%[0 1001];
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
    ...& numMatchedTargets >= 3 ...
    ...& ensembleOneSecond ... %cuts off a lot of the earlier
    & numCellsEachEns==1 ...
    ...& ensDate >= 211101 ...
    ...& outVars.hzEachEns == 10 ...
    ...& outVars.hzEachEns >= 9 & outVars.hzEachEns <= 12 ...
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


%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials ...
        & All(ind).out.exp.lowRunTrials ...
        & All(ind).out.exp.stimSuccessTrial ...
        & (All(ind).out.exp.visID == 1 |  All(ind).out.exp.visID == 0 ) ...
            ;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end

% %% Orientation Tuning and OSI
% 
% [All, outVars] = getTuningCurve(All, opts, outVars);
% [All, outVars] = calcOSI(All, outVars);
% [All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
% [All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically

%% Distance of Ensemble
[All, outVars] = defineDistanceTypes(All, outVars);
opts.distType = 'min';
[outVars] = grandDistanceMaker(opts,All,outVars);


%% Plot Mean Distance Responses

plotEnsembleDistanceResponse(outVars,100,1)

%% Plot Distance Plots

opts.distBins = 0:25:1000; %must be set to match popDist
plotResponseByDistance(outVars,opts);

%% Compare Distance responses
figure(102);clf

dataInPlots=[];
distTypes = {'min' 'geo' 'mean' 'harm' 'median' 'centroid'};
for i =[1 3]; %1:6
    disp(['working on ' distTypes{i}])
    opts.distType = distTypes{i}; %options: min geo mean harm
    opts.distBins = 0:10:350; %can be set variably 0:25:1000 is defaultt
    CellToUseVar = 'anal.cellsToInclude';
    [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
    ax = subplot(2,3,i);
    opts.distAxisRange = [0 350]; %[0 350] is stand
    [eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
    dataInPlots{i}=outDat{1};
    eHandle{1}.CapSize =0;
    title(distTypes{i})
    drawnow
end
disp('done')

%% Just a few with different binning
figure(103);clf;
dataInPlots =[];

ax = subplot(1,2,1);
opts.distType = 'min';
opts.distBins = 0:10:350; %can be set variably 0:25:1000 is defaultt
opts.distAxisRange = [0 250]; %[0 350] is stand
CellToUseVar = 'anal.cellsToInclude';
[popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
[eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
dataInPlots{1}=outDat{1};
eHandle{1}.CapSize =0;
title('min')

ax = subplot(1,2,2);
opts.distType = 'mean';
opts.distBins = 0:10:500; %can be set variably 0:25:1000 is defaultt
opts.distAxisRange = [0 450]; %[0 350] is stand
CellToUseVar = 'anal.cellsToInclude';
[popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
[eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
dataInPlots{2}=outDat{1};
eHandle{1}.CapSize =0;
title('mean')
ylim([-0.075 0.075])


%% create new


for ind = 1:numExps
    if ~ismember(ind,noOpsinDataInd)
        vals = All(ind).out.info.value1020;
        defPosThresh = prctile(vals,70);
        All(ind).out.info.defPos = vals>=defPosThresh;
        
        defNegThresh = prctile(vals,30);
        All(ind).out.info.defNeg = vals<=defNegThresh;
    else
        All(ind).out.info.defPos = zeros([1 All(ind).out.anal.numCells]);
        All(ind).out.info.defNeg = zeros([1 All(ind).out.anal.numCells]);
        All(ind).out.info.opsinNegative = zeros([1 All(ind).out.anal.numCells]);
    end
end

backupEnsemblesToUse = outVars.ensemblesToUse;

%% to match  other plots
figure(104);clf
hold on
ax = subplot(1,1,1);

opts.distBins = 0:5:50; %00:40:350; %0:25:350; %can be set variably 0:25:1000 is defaultt
opts.distType = 'min';
opts.distAxisRange = [0 250]; %[0 350] is stand

% backupEnsemblesToUse = outVars.ensemblesToUse;
% noUnstimableCount = find(countUSC==0);
%  limEnsembleToUse = outVars.ensemblesToUse & ~ismember(outVars.ensIndNumber,[1 4]);
%  outVars.ensemblesToUse = limEnsembleToUse;
disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
% 
CellToUseVar =[];
[popRespDistDefault] = popDistMaker(opts,All,CellToUseVar,0);
p1 = plotDistRespGeneric(popRespDistDefault,outVars,opts,ax);
p1{1}.Color=rgb('black');
p1{1}.CapSize = 0;
outVars.ensemblesToUse = backupEnsemblesToUse;
hold on
drawnow
ylim([-0.05 0.11])


% backupEnsemblesToUse = outVars.ensemblesToUse;
% noUnstimableCount = find(countUSC==0);
%  limEnsembleToUse = outVars.ensemblesToUse & ~ismember(outVars.ensIndNumber,[1 4]);
%  outVars.ensemblesToUse = limEnsembleToUse;
disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
% 
CellToUseVar ='info.defPos';%[];
[popRespDistDefault] = popDistMaker(opts,All,CellToUseVar,0);
p1 = plotDistRespGeneric(popRespDistDefault,outVars,opts,ax);
p1{1}.Color=rgb('red');
p1{1}.CapSize = 0;
outVars.ensemblesToUse = backupEnsemblesToUse;
hold on
drawnow
ylim([-0.05 0.11])

% backupEnsemblesToUse = outVars.ensemblesToUse;
% noUnstimableCount = find(countUSC==0);
%  limEnsembleToUse = outVars.ensemblesToUse & ismember(outVars.ensIndNumber,[1 2 5]);
%  outVars.ensemblesToUse = limEnsembleToUse;
% disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
% 
CellToUseVar ='info.defNeg'; %'info.opsinNegative';% [];
[popRespDistData] = popDistMaker(opts,All,CellToUseVar,0);
p1 = plotDistRespGeneric(popRespDistData,outVars,opts,ax);
p1{1}.Color=rgb('ForestGreen');
p1{1}.CapSize = 0;
outVars.ensemblesToUse = backupEnsemblesToUse;
hold on
drawnow
ylim([-0.05 0.25])

%% Tick through each experiment

figure(104);clf
hold on
ax = subplot(1,1,1);

opts.distBins = 00:10:350; %0:25:350; %can be set variably 0:25:1000 is defaultt
opts.distType = 'min';
opts.distAxisRange = [0 250]; %[0 350] is stand

outVars.ensemblesToUse = backupEnsemblesToUse;

% backupEnsemblesToUse = outVars.ensemblesToUse;
% noUnstimableCount = find(countUSC==0);
%  limEnsembleToUse = outVars.ensemblesToUse & ~ismember(outVars.ensIndNumber,[1 4]);
%  outVars.ensemblesToUse = limEnsembleToUse;
disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
% 
CellToUseVar =[];
[popRespDistDefault] = popDistMaker(opts,All,CellToUseVar,0);
p1 = plotDistRespGeneric(popRespDistDefault,outVars,opts,ax);
p1{1}.Color=rgb('black');
p1{1}.CapSize = 0;
outVars.ensemblesToUse = backupEnsemblesToUse;
hold on
drawnow
ylim([-0.05 0.11])

colors = colorMapPicker(numel(All),'viridis');
for i=1:max(outVars.ensIndNumber); %5
%     backupEnsemblesToUse = outVars.ensemblesToUse;
    % noUnstimableCount = find(countUSC==0);
    limEnsembleToUse = outVars.ensemblesToUse & ismember(outVars.ensIndNumber,[i]);
    if sum(limEnsembleToUse)>0
        outVars.ensemblesToUse = limEnsembleToUse;
        disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
        %
        CellToUseVar ='info.defNeg'; %info.opsinNegative';% [];
        [popRespDistDefault] = popDistMaker(opts,All,CellToUseVar,0);
        p1 = plotDistRespGeneric(popRespDistDefault,outVars,opts,ax);
        p1{1}.Color=colors{i}; %rgb('forestgreen');
        p1{1}.CapSize = 0;
       
        hold on
        ylim([-0.05 0.25])
        drawnow
        disp('press the any key')
        pause
    end
     outVars.ensemblesToUse = backupEnsemblesToUse;
end


%%
%% ID opsin negative offline
figure(104);clf
hold on
ax = subplot(1,1,1);

opts.distBins = 0:5:50; 00:10:350; %0:25:350; %can be set variably 0:25:1000 is defaultt
opts.distType = 'min';
opts.distAxisRange = [0 250]; %[0 350] is stand

% backupEnsemblesToUse = outVars.ensemblesToUse;
% noUnstimableCount = find(countUSC==0);
%  limEnsembleToUse = outVars.ensemblesToUse & ~ismember(outVars.ensIndNumber,[1 4]);
%  outVars.ensemblesToUse = limEnsembleToUse;
disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
% 
CellToUseVar =[];
[popRespDistDefault] = popDistMaker(opts,All,CellToUseVar,0);
p1 = plotDistRespGeneric(popRespDistDefault,outVars,opts,ax);
p1{1}.Color=rgb('black');
p1{1}.CapSize = 0;
outVars.ensemblesToUse = backupEnsemblesToUse;
hold on
drawnow
ylim([-0.05 0.11])




% backupEnsemblesToUse = outVars.ensemblesToUse;
% noUnstimableCount = find(countUSC==0);
 limEnsembleToUse = outVars.ensemblesToUse;% & ismember(outVars.ensIndNumber,[3]);
 outVars.ensemblesToUse = limEnsembleToUse;
disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
% 
CellToUseVar ='info.opsinNegative';% [];
[popRespDistData] = popDistMaker(opts,All,CellToUseVar,0);
p1 = plotDistRespGeneric(popRespDistData,outVars,opts,ax);
p1{1}.Color=rgb('ForestGreen');
p1{1}.CapSize = 0;
outVars.ensemblesToUse = backupEnsemblesToUse;
hold on
drawnow
ylim([-0.05 0.11])

%% Compare Opsin Positive between different Expression Types

ExprUsed = outVars.indExpressionType(IndsUsed);
outVars.uniqueExpressionTypes(ExprUsed)

