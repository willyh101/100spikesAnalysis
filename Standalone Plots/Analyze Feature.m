%% Load Experiments/Setup
clear
close all

masterTic = tic;

addpath(genpath('100spikesAnalysis'))
%% loadLists

oriLoadList;
% oriLoadListOnlyInt;

% % allLoadList;

% loadPath = 'path/to/outfiles/directory';
% loadPath = 'T:\Outfiles';

loadPath = 'C:\Users\ian\Dropbox\Outfiles';

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

%% Rematch Targeted Cells
opts.matchDistanceMicrons = 12; 15.6; %6;
rematchTargetedCells;

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
opts.minMeanThreshold = -0.25;
opts.maxMeanThreshold = inf;

opts.verbose =0;
[All, cellExcludeResults] = cellExcluder(All,opts); 
allResults = cat(1,cellExcludeResults{:});
disp(['In total ' num2str(sum(allResults)) ' Cells Excluded. ' num2str(mean(allResults)*100,2) '%']);
disp(['Overall ' num2str(sum(~allResults)) ' Cells Passed!'])

opts.minNumCellsInd= 250;
tooFewCellsInds = cellfun(@(x) sum(~x)<opts.minNumCellsInd,cellExcludeResults);
disp([ num2str(sum(tooFewCellsInds)) ' inds have < ' num2str(opts.minNumCellsInd) ' cells, and should be exccluded']);

%% Make all dataPlots into matrixes of mean responses
%%Determine Vis Responsive and Process Correlation

opts.visAlpha = 0.05;
muPerPx=800/512;
%oftarget risk params
opts.thisPlaneTolerance = 15/muPerPx;11.25;%7.5;%1FWHM%10; %in pixels
opts.onePlaneTolerance = 30/muPerPx;22.5;%15;%2FWHM %20;
opts.distBins =  [0:25:1000]; [0:25:1000];
opts.skipVis =1;

[All, outVars] = meanMatrixVisandCorr(All,opts,outVars); %one of the main analysis functions

visPercent = outVars.visPercent;
outVars.visPercentFromExp = visPercent;
ensIndNumber =outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber; %in order of uniqueStims


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

% numTrialsPerEns =[];numTrialsPerEnsTotal=[]; numTrialsNoStimEns=[];
% for ind=1:numExps
%     us=unique(All(ind).out.exp.stimID);
% 
%     for i=1:numel(us)
%         trialsToUse = All(ind).out.exp.lowMotionTrials &...
%             All(ind).out.exp.lowRunTrials &...
%             All(ind).out.exp.stimSuccessTrial &...
%             All(ind).out.exp.stimID == us(i) & ...
%             (All(ind).out.exp.visID == 1 | All(ind).out.exp.visID == 0); %restrict just to no vis stim conditions
% 
%         numTrialsPerEns(end+1)=sum(trialsToUse);
%         numTrialsPerEnsTotal(end+1) = sum(All(ind).out.exp.stimID == us(i));
% 
%         if i==1
%             numTrialsNoStimEns(ind) = sum(trialsToUse);
%         end
%     end
% end
% numTrialsPerEns(numSpikesEachStim==0)=[];
% numTrialsPerEnsTotal(numSpikesEachStim==0)=[];
% 
% %ID inds to be excluded
% opts.IndsVisThreshold = 0.1;0.05; %default 0.05
% 
% highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<opts.IndsVisThreshold)); %remove low vis responsive experiments
% lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));
% lowCellCount = ismember(ensIndNumber,find(tooFewCellsInds));
% 
% 
% %exclude certain expression types:
% uniqueExpressionTypes = outVars.uniqueExpressionTypes;
% excludedTypes ={'AAV CamK2' 'Ai203' 'neo-IV Tre 2s' 'IUE CAG' 'SepW1 CAG 2s' };
% 
% 
% exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
% excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);
% 
% % only include times where rate == numpulses aka the stim period is 1s.
% ensembleOneSecond = outVars.numSpikesEachEns./outVars.numCellsEachEns == outVars.hzEachEns;
% 
% %spot to add additional Exclusions
% excludeInds = ismember(ensIndNumber,[]); %Its possible that the visStimIDs got messed up
% % excludeInds = ismember(ensIndNumber,[]); 
% 
% %Options
% opts.numSpikeToUseRange = [90 110];[1 inf];[80 120];%[0 1001];
% opts.ensStimScoreThreshold = 0.5; % default 0.5
% opts.numTrialsPerEnsThreshold = 5; % changed from 10 by wh 4/23 for testing stuff
% 
% lowBaseLineTrialCount = ismember(ensIndNumber,find(numTrialsNoStimEns<opts.numTrialsPerEnsThreshold));
% 
% 
% ensemblesToUse = numSpikesEachEns > opts.numSpikeToUseRange(1) ...
%     & numSpikesEachEns < opts.numSpikeToUseRange(2) ...
%     & highVisPercentInd ...
%     & lowRunInds ...
%     & ensStimScore > opts.ensStimScoreThreshold ... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
%     & ~excludeInds ...
%     & numTrialsPerEns > opts.numTrialsPerEnsThreshold ... ;%10;%&...
%     & ~lowBaseLineTrialCount ...
%     & ~ensHasRed ...
%     & ~excludeExpressionType ...
%     & ~ensMissedTarget ...
%     & numMatchedTargets >= 3 ...
%     & ensembleOneSecond ... %cuts off a lot of the earlier
%     & numCellsEachEns==10 ...
%     ... & ensDate > 210420 ...
%     ...& outVars.hzEachEns == 10 ...
%     ...& outVars.hzEachEns >= 9 & outVars.hzEachEns <= 12 ...
%     & ~lowCellCount ...
%     ;
% 
% %%remove repeats
%  [ensemblesToUse, outVars] = removeRepeatsFromEnsemblesToUse(ensemblesToUse,outVars);
% 
% indsSub = ensIndNumber(ensemblesToUse);
% IndsUsed = unique(ensIndNumber(ensemblesToUse));
% 
% sum(ensemblesToUse)
% 
% outVars.ensemblesToUse      = ensemblesToUse;
% outVars.IndsUsed            = IndsUsed;
% outVars.indsSub             = indsSub;
% outVars.numTrialsPerEns     = numTrialsPerEns;
% outVars.highVisPercentInd    = highVisPercentInd;
% outVars.lowRunInds           = lowRunInds;
% 
% %%Optional: Where are the losses comming from
% 
% disp(['Fraction of Ens correct Size: ' num2str(mean(numSpikesEachEns > opts.numSpikeToUseRange(1) & numSpikesEachEns < opts.numSpikeToUseRange(2)))]);
% disp(['Fraction of Ens highVis: ' num2str(mean(highVisPercentInd))]);
% disp(['Fraction of Ens lowRun: ' num2str(mean(lowRunInds))]);
% disp(['Fraction of Ens high stimScore: ' num2str(mean(ensStimScore>opts.ensStimScoreThreshold))]);
% disp(['Fraction of Ens high trial count: ' num2str(mean(numTrialsPerEns>opts.numTrialsPerEnsThreshold))]);
% disp(['Fraction of Control Ens high trial count: ' num2str(mean(~lowBaseLineTrialCount))]);
% disp(['Fraction of Ens No ''red'' cells shot: ' num2str(mean(~ensHasRed))]);
% disp(['Fraction of Ens usable Expression Type: ' num2str(mean(~excludeExpressionType))]);
% disp(['Fraction of Ens enough targets detected by s2p: ' num2str(mean(~ensMissedTarget))]);
% disp(['Fraction of Ens number targets matched >=3: ' num2str(mean(numMatchedTargets >= 3))]);
% disp(['Fraction of Ens Stim took 1s (aka correct stim Rate): ' num2str(mean(ensembleOneSecond))]);
% disp(['Fraction of Ens that were not repeats: ' num2str(mean(~outVars.removedRepeats)) ]);
% disp(['Fraction of Ens high Cell Count: ' num2str(mean(~lowCellCount))]);
% 
% 
% disp(['Total Fraction of Ens Used: ' num2str(mean(ensemblesToUse))]);
% % disp([num2str(sum(ensemblesToUse)) ' Ensembles Included'])
% disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
% disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
% disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);

StandardEnsemblesToUse

%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials ...
        & All(ind).out.exp.lowRunTrials ...
        & All(ind).out.exp.stimSuccessTrial ...
        & (All(ind).out.exp.visID == 1 |  All(ind).out.exp.visID == 0 ) ...
            ;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end

%% Orientation Tuning and OSI

[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
[All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically

%% Distance of Ensemble
[All, outVars] = defineDistanceTypes(All, outVars);
opts.distType = 'min';
[outVars] = grandDistanceMaker(opts,All,outVars);

%% Plot Mean Responses
opts.ensOSImethod = 'ensOSI';
opts.plotFit=1;
opts.figToUse1 = 44;
opts.defaultColor = rgb('firebrick');
result = plotPopResponseEnsOSI(outVars, opts);

legend off
disp([opts.ensOSImethod '. Fit R2= ' num2str(result.Rsq,3) ' pVal= ' num2str(result.pVal(1),3) ])
ylim([-0.15 0.05])
yticks(-.15:.05:.05)


%%
opts.ensOSImethod = 'meanEnsOSI';
opts.plotFit=1;
opts.figToUse1 = 45;
opts.defaultColor = rgb('firebrick');

result = plotPopResponseEnsOSI(outVars, opts);
legend off
disp([opts.ensOSImethod '. Fit R2= ' num2str(result.Rsq,3) ' pVal= ' num2str(result.pVal(1),3) ])
ylim([-0.15 0.05])
yticks(-.15:.05:.05)

%%
outVars.EnsOSIxMeanOSI = outVars.ensOSI .* outVars.meanEnsOSI;
    
opts.ensOSImethod = 'EnsOSIxMeanOSI';
opts.plotFit=1;
opts.figToUse1 = 46;
opts.defaultColor = rgb('firebrick');

result = plotPopResponseEnsOSI(outVars, opts);
legend off
disp([opts.ensOSImethod '. Fit R2= ' num2str(result.Rsq,3) ' pVal= ' num2str(result.pVal(1),3) ])
ylim([-0.15 0.05])
yticks(-.15:.05:.05)


%  %%
% outVars.EnsMaxDxMeanOSI = outVars.ensMaxD .* outVars.ensOSI;
%     
% opts.ensOSImethod = 'ensMaxD'; %'EnsMaxDxMeanOSI';
% opts.plotFit=1;
% opts.figToUse1 = 47;
% opts.defaultColor = rgb('firebrick');
% 
% result = plotPopResponseEnsOSI(outVars, opts);
% legend off
% disp([opts.ensOSImethod '. Fit R2= ' num2str(result.Rsq,3) ' pVal= ' num2str(result.pVal(1),3) ])
% ylim([-0.15 0.05])

%% Plot Population Response read out in different cells


opts.ensemblesToPlot = outVars.ensemblesToUse; 
opts.useVisCells =0;
opts.useTunedCells =0;
opts.minNumberOfCellsPerCondition = -1;

%unTuned
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<0.3 &  outVars.meanEnsOSI<0.5;
sum(opts.ensemblesToPlot)
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05';
[UnTunedVis] = subsetPopResponse(All,outVars,opts);
opts.variableCellFun =  'outVars.pVisR{ind} > 0.1';
[UnTunedNotVis] = subsetPopResponse(All,outVars,opts);

opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<0.3 &  outVars.meanEnsOSI>0.5; 
sum(opts.ensemblesToPlot)

opts.variableCellFun =  'outVars.pVisR{ind} < 0.05';
[MxTunedVis] = subsetPopResponse(All,outVars,opts);
opts.variableCellFun =  'outVars.pVisR{ind} > 0.1';
[MxTunedNotVis] = subsetPopResponse(All,outVars,opts);

opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI>0.7 &  outVars.meanEnsOSI>0.5; 
sum(opts.ensemblesToPlot)

opts.variableCellFun =  'outVars.pVisR{ind} < 0.05';
[CoTunedVis] = subsetPopResponse(All,outVars,opts);
opts.variableCellFun =  'outVars.pVisR{ind} > 0.1';
[CoTunedNotVis] = subsetPopResponse(All,outVars,opts);


figure(48);clf;hold on
dat = {UnTunedNotVis, UnTunedVis, [], MxTunedNotVis, MxTunedVis, [], CoTunedNotVis, CoTunedVis};
distColors = {rgb('grey') rgb('purple'), rgb('black'), rgb('grey') rgb('purple'), rgb('black'), rgb('grey') rgb('purple')};
p = plotSpread(dat,'distributionColors',distColors,'ShowMM',4);
r = refline(0);
r.Color= rgb('grey');
r.LineStyle = ':';

xticks([1.5 4.5 7.5])
xticklabels({'Untuned' 'Mixed Tuned' 'Co Tuned'})
ylabel('Mean Evoked dF/F')

%% Plot Mean Population Response in different Ori Bands;

oriVals = [NaN 0:45:315];
diffsPossible = [0 45 90 135 180];
plotOrientation=1;

for i = 1:numel(outVars.ensemblesToUse);
    ind = ensIndNumber(i);
    
    cellOris = oriVals(outVars.prefOris{ind});
    cellOrisDiff = abs(cellOris-outVars.ensPO(i));
    cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
    %
    if plotOrientation
        cellOrisDiff(cellOrisDiff==135)=45;
        cellOrisDiff(cellOrisDiff==180)=0;
        diffsPossible = [0 45 90];
    end
    
    outVars.ensOriDiff{i} = cellOrisDiff;
end


opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI>0.7 &  outVars.meanEnsOSI>0.5; 
sum(opts.ensemblesToPlot)
opts.useVisCells =1;
opts.useTunedCells =1;
opts.minNumberOfCellsPerCondition =-20; %set to -1 to ignore


opts.variableCellFun =  'outVars.ensOriDiff{i} == 0';
[OriDiff0] = subsetPopResponse(All,outVars,opts);
opts.variableCellFun =  'outVars.ensOriDiff{i} == 45';
[OriDiff45] = subsetPopResponse(All,outVars,opts);
opts.variableCellFun =  'outVars.ensOriDiff{i} == 90';
[OriDiff90] = subsetPopResponse(All,outVars,opts);
opts.variableCellFun =  'outVars.ensOriDiff{i} == 135';
[OriDiff135] = subsetPopResponse(All,outVars,opts);
opts.variableCellFun =  'outVars.ensOriDiff{i} == 180';
[OriDiff180] = subsetPopResponse(All,outVars,opts);


figure(48);clf;hold on
dat = {OriDiff0,OriDiff45,OriDiff90,OriDiff135,OriDiff180};
distColors = colorMapPicker(5,'plasma');

p = plotSpread(dat,'distributionColors',distColors,'ShowMM',4);
r = refline(0);
r.Color= rgb('grey');
r.LineStyle = ':';
% 
% xticks([1.5 4.5 7.5])
% xticklabels({'Untuned' 'Mixed Tuned' 'Co Tuned'})
% ylabel('Mean Evoked dF/F')

if plotOrientation==1;
xlim([0.5 3.5])

xticklabels({['Iso (0' char(176) ')'] [char(177) '45' char(176)] ['Ortho (90' char(176) ')'] })

else
    xticklabels({['Iso (0' char(176) ')']...
        [char(177) '45' char(176)]...
        ['Ortho (' char(177) '90' char(176) ')']...
        [char(177) '135' char(176)]...
        ['Reverse (180' char(176) ')']...
        })
end


ylabel('Mean Evoked dF/F')

ranksum(OriDiff0,OriDiff90)
%% Plot Split by Mean OSI
opts.distType = 'min';
opts.distBins =[15:15:150];[0:25:150]; 
opts.distAxisRange = [0 150];
opts.plotTraces = 0; 
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<0.3;
opts.criteriaToSplit = outVars.meanEnsOSI;
opts.criteriaBins = [0 0.4 0.5 inf];
% opts.useVisAndTunedCells =1; 
opts.useVisCells =0;
opts.useTunedCells =0; %don't use tuned without vis

plotDistByCriteria(All,outVars,opts,15)
figure(16);
ylim([-0.15 0.15]);

%%  Plot Split by Ensemble OSI
opts.distType = 'min';
opts.distBins =[0:25:150]; 
opts.plotTraces = 0; 
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.meanEnsOSI>0.5;
opts.criteriaToSplit = outVars.ensOSI;
opts.criteriaBins = [0 0.3 0.7 inf];
% opts.useVisAndTunedCells =1; 
opts.useVisCells =1;
opts.useTunedCells =0; %don't use tuned without vis

plotDistByCriteria(All,outVars,opts,17)
figure(18);
ylim([-0.15 0.1]);

%% Plot Dist by Two Criteria
opts.distType = 'min';
opts.distBins =[0:25:150]; 
opts.plotTraces = 0; 
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.criteriaName = 'ensOSI';
opts.criteria = outVars.ensOSI;
opts.criteriaBins = [0 0.3 0.7 inf];

opts.criteria2Name = 'meanOSI';
opts.criteria2 = outVars.meanEnsOSI;
opts.criteria2Bins = [0 0.5 inf];
% opts.useVisAndTunedCells =1; 
opts.useVisCells =1;
opts.useTunedCells =0; %don't use tuned without vis

plotDistByTwoCriteria(All,outVars,opts,13)
figure(13)
ylim([-0.15 0.1]);



%%
muPerPx = 800/512;
opts.thisPlaneTolerance =15/muPerPx;% 15/muPerPx;
opts.onePlaneTolerance = 30/muPerPx; %30/muPerPx;

recalcOffTargetRisk;
%% Plot Distance Curves in different Ori Bands

opts.distType = 'min';
opts.distBins =15:15:150; 0:25:150; %15:15:150; %[0:25:150];%[0:200:1000];%[0:25:150]; [0:100:1000];% [0:25:400];


opts.plotOrientation =1;%as opposed to Direction
opts.minNumberOfCellsPerCondition =-10; %set to -1 to ignore
opts.ensemblesToPlot = outVars.ensemblesToUse  & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI>0.7 &  outVars.meanEnsOSI>0.5;
% opts.ensemblesToPlot = outVars.ensemblesToUse  & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI>0.7 &  outVars.meanEnsOSI>0.5 & outVars.ensMaxD>500;
% opts.ensemblesToPlot = outVars.ensemblesToUse  & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI>0.7 &  outVars.meanEnsOSI>0.5 & outVars.ensMaxD<400;

% opts.ensemblesToPlot = outVars.ensemblesToUse  & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<0.3 &  outVars.meanEnsOSI<0.5;

 plotDistByOri(All,outVars,opts)
% plotDistByClosestTargetOri(All,outVars,opts)

 figure(10);
ylim([-0.1 0.125])


figure(14);
ylim([-0.225 0.5])


%%Plot Space and Feature V3 
opts.distType = 'min';
% opts.distBins =15:15:150; %15:15:150; %[15:15:150]; [opts.thisPlaneTolerance*muPerPx:10:150]; %[0:25:150]; 
opts.plotTraces = 0;
opts.useVisCells = 0;
opts.useTunedCells = 0; %don't use tuned without vis
figure(120); clf
lim =[-0.1 0.125]; [-0.4 0.25];

highMeanThresh = 0.5; %0.5; % 0.4687;
lowMeanThresh = 0.5;

lowEnsThresh = 0.3;
highEnsThresh = 0.7;

closeVal =  400
farVal =500

outInfo=[];
axs = [];
ax = subplot(3,2,1);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<lowEnsThresh & outVars.meanEnsOSI<lowMeanThresh;
opts.criteriaToSplit =  outVars.ensMaxD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('grey');
axs = [axs ax];

ax = subplot(3,2,2);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('grey');
axs = [axs ax];


ax = subplot(3,2,3);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<lowEnsThresh & outVars.meanEnsOSI>highMeanThresh;
% opts.criteriaToSplit = outVars.ensMaxD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('sienna');
axs = [axs ax];

ax = subplot(3,2,4);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('sienna');
axs = [axs ax];


ax = subplot(3,2,5);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI>highEnsThresh & outVars.meanEnsOSI>highMeanThresh;
% opts.criteriaToSplit = outVars.ensMaxD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('Magenta');
axs = [axs ax];

ax = subplot(3,2,6);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('Magenta');
axs = [axs ax];

linkaxes(axs)
ylim(lim);

disp('pVal first point diff from zero')
for i =1:6
    disp(num2str(signrank(outInfo{i}{1}.dat(:,1),0)))
end

[p h] = ranksum(outInfo{5}{1}.dat(:,1),outInfo{6}{1}.dat(:,1));
disp(['Tuned Near vs Far p= ' num2str(p)]);

[p h] = ranksum(outInfo{1}{1}.dat(:,1),outInfo{5}{1}.dat(:,1));
disp(['Near Untuned vs Tuned p= ' num2str(p)]);

[p h] = ranksum(outInfo{2}{1}.dat(:,1),outInfo{6}{1}.dat(:,1));
disp(['Far Untuned vs Tuned p= ' num2str(p)]);

%%
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<0.7;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;
opts.variableCellFun =  '(outVars.pVisR{ind} < 0.05 & outVars.distToEnsemble{i}<25)';
[closeExcitation] = subsetPopResponse(All,outVars,opts);

figure(19);clf
scatter(outVars.meanEnsOSI,closeExcitation,[],rgb('sienna'),'filled')
refline(0)
xlabel('MeanEnsOsi') ;

x = outVars.meanEnsOSI(opts.ensemblesToPlot)';
y = closeExcitation(opts.ensemblesToPlot);
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp(['Pval is: ' num2str(pVal(1))])

%%
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.meanEnsOSI>0.5;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;
opts.variableCellFun =  '(outVars.pVisR{ind} < 0.05 & outVars.distToEnsemble{i}<30)';
[closeExcitation] = subsetPopResponse(All,outVars,opts);

figure(19);clf
scatter(outVars.ensOSI,closeExcitation,[],rgb('sienna'),'filled')
refline(0)
xlabel('EnsOsi') ;

x = outVars.ensOSI(opts.ensemblesToPlot)';
y = closeExcitation(opts.ensemblesToPlot);
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp(['Pval is: ' num2str(pVal(1))])

%% cut older analyses
%


%% Plot Space and Feature
% 
% opts.distType = 'min';
% opts.distBins =[0:25:150];%[15:20:150];% [0:25:400];
% opts.distAxisRange = [min(opts.distBins) max(opts.distBins)]; %[0 350];
% opts.plotTraces =0;
% opts.useVisCells =1;
% opts.useTunedCells =0; %don't use tuned without vis
% 
% plotSpaceAndFeature(All,outVars,opts,11)
% ylim([-0.15 0.15]);
% 
% 
% %% Plot Space and  Feature v2
% %% Plot Dist by Two Criteria
% opts.distType = 'min';
% opts.distBins =[15:15:150]; [0:25:150]; 
% opts.plotTraces = 0; 
% opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
% opts.criteria2Name = 'ensOSI';
% opts.criteria2 = outVars.ensOSI;
% opts.criteria2Bins = [0 0.3 0.7 inf]
% 
% opts.criteriaName = 'Spread';
% opts.criteria = outVars.ensMaxD;
% opts.criteriaBins = [0 400 500 inf];
% opts.useVisCells =1; 
% opts.useTunedCells =0; 
% 
% plotDistByTwoCriteria(All,outVars,opts,21)
% % figure(13)
% % ylim([-0.15 0.1]);




%%
toc(masterTic)