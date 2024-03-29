%% Load Experiments/Setup
clear
close all

masterTic = tic;

addpath(genpath('100spikesAnalysis'))
%% loadLists

oriLoadList_relevant;
GMNLoadList;
% % allLoadList;

% loadPath = 'path/to/outfiles/directory';
% loadPath = 'T:\Outfiles';
% loadPath = '/Users/gregoryhandy/Research_Local/outputdata1';
% loadPath = '/Users/greghandy/Research_Local_v2/'; % where ever your files are
loadPath = '/Users/ianoldenburg/Dropbox/Outfiles_Pulled_221209';

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
muPerPx=800/512;

%oftarget risk params
opts.thisPlaneTolerance = 15/muPerPx;11.25;%7.5;%1FWHM%10; %in pixels
opts.onePlaneTolerance = 30/muPerPx;22.5;%15;%2FWHM %20;
opts.distBins =  [0:25:1000]; [0:25:1000];

opts.skipVis =0;

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
opts.IndsVisThreshold = 0.05;0.05; %default 0.05

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<opts.IndsVisThreshold)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));
lowCellCount = ismember(ensIndNumber,find(tooFewCellsInds));


%exclude certain expression types:
uniqueExpressionTypes = outVars.uniqueExpressionTypes;
excludedTypes ={'AAV CamK2' 'Ai203' 'neo-IV Tre 2s' 'IUE CAG' 'PHP Tre' 'control' };%'SepW1 CAG 2s'


exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

% only include times where rate == numpulses aka the stim period is 1s.
ensembleOneSecond = outVars.numSpikesEachEns./outVars.numCellsEachEns == outVars.hzEachEns;



%how many ensemblesPer Ind
uEIN =unique(outVars.ensIndNumber);
indEnsCount =[];
for i=1:numel(uEIN)
    e = uEIN(i);
    indEnsCount(i) = sum(outVars.ensIndNumber==e);
end

%spot to add additional Exclusions
%  excludeInds = ismember(ensIndNumber,find(indEnsCount>50)); %Its possible that the visStimIDs got messed up
excludeInds = ismember(ensIndNumber,[]); 

%Options
opts.numSpikeToUseRange = [90 110];[1 inf];[80 120];%[0 1001];
opts.ensStimScoreThreshold = 0.5; % default 0.5
opts.numTrialsPerEnsThreshold = 10; %5 changed from 10 by wh 4/23 for testing stuff

lowBaseLineTrialCount = ismember(ensIndNumber,find(numTrialsNoStimEns<opts.numTrialsPerEnsThreshold));


ensemblesToUse = numSpikesEachEns > opts.numSpikeToUseRange(1) ...
    & numSpikesEachEns < opts.numSpikeToUseRange(2) ...
    ...& highVisPercentInd ...
    & lowRunInds ...
    & ensStimScore > opts.ensStimScoreThreshold ... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    & ~excludeInds ...
    ...& numTrialsPerEns > opts.numTrialsPerEnsThreshold ... ;%10;%&...
    ...& ~lowBaseLineTrialCount ...
    ...& ~ensHasRed ...
    & ~excludeExpressionType ...
    & ~ensMissedTarget ...
    & numMatchedTargets >= 7 ...
    ...& ensembleOneSecond ... %cuts off a lot of the earlier
    & numCellsEachEns==10 ...
    ...& (ensDate < 220000 & ensDate>210420) ...
    & ensDate < 220000 ...
    ...& ensDate < 210000 ...
    ...& (ensDate > 220000 & ensDate>220226) ...
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
disp(['Fraction of Ens number targets matched >=7: ' num2str(mean(numMatchedTargets >= 3))]);
disp(['Fraction of Ens Stim took 1s (aka correct stim Rate): ' num2str(mean(ensembleOneSecond))]);
disp(['Fraction of Ens that were not repeats: ' num2str(mean(~outVars.removedRepeats)) ]);
disp(['Fraction of Ens high Cell Count: ' num2str(mean(~lowCellCount))]);


disp(['Total Fraction of Ens Used: ' num2str(mean(ensemblesToUse))]);
% disp([num2str(sum(ensemblesToUse)) ' Ensembles Included'])
disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);

% StandardEnsemblesToUse

%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials ...
        & All(ind).out.exp.lowRunTrials ...
        & All(ind).out.exp.stimSuccessTrial ...
        & (All(ind).out.exp.visID == 1 |  All(ind).out.exp.visID == 0 ) ...
            ;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end

%% basic funs
outVars.defaultColorMap = 'viridis';
% plotAllEnsResponse(outVars)
plotResponseBySize(outVars,0)
% plotPopResponseBySession(All,outVars)
% plotPopResponseByExpressionType(All,outVars);
% [All, outVars] = createTSPlotByEnsSize(All,outVars);
% [All, outVars] = createTSPlotByEnsSizeAllVis(All,outVars);


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


%% Plot Distance 
figure(1003);clf
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
eHandle{1}.CapSize =0;
title('min')

%% Plot Contrast responses
opts.ensemblesToPlot = outVars.ensemblesToUse; 
opts.useVisCells =0;
opts.useTunedCells =0;
opts.minNumberOfCellsPerCondition = -1;
% opts.variableCellFun =  'outVars.pVisR{ind} < 0.05';

mCRFnoStim = []; mCRFens = [];
eCRFnoStim = []; eCRFens = [];
for i =1:6
    opts.visCond=i;

[test1] = subsetPopResponseMatchedNoStim(All,outVars,opts);
[test2] = subsetPopResponse(All,outVars,opts);

mCRFnoStim(i) = mean(test1(opts.ensemblesToPlot));
mCRFens(i) = mean(test2(opts.ensemblesToPlot));

eCRFnoStim(i) = ste(test1(opts.ensemblesToPlot));
eCRFens(i) = ste(test2(opts.ensemblesToPlot));
end
opts.visCond=1;

figure(707);clf
hold on
e1= errorbar(mCRFnoStim,eCRFnoStim);
e1.Color = rgb('grey');
e1.LineStyle=':';
e1.LineWidth = 2;

e2 = errorbar(mCRFens,eCRFens);
e2.Color = rgb('FireBrick');
% e1.LineStyle=':';
e2.LineWidth = 2;

xticks(1:6)
xticklabels({'0%', '1%', '4%', '10%', '40%', '100%'})
xlabel('Contrast')
ylabel('Pop Response (Evoked \DeltaF/F)')


%% Dist Response during Contrasts
figure(1003);clf
hold on
dataInPlots =[];

ax = subplot(1,1,1);
opts.distType = 'min';
opts.distBins = 0:10:250;%15:15:250; %can be set variably 0:25:1000 is defaultt
opts.distAxisRange = [0 250]; %[0 350] is stand
% CellToUseVar = 'anal.cellsToInclude & All(ind).out.anal.pVisR<0.05';%[];
CellToUseVar = ''; %'anal.cellsToInclude ';%[];

clist = colorMapPicker(6,'plasma');
clist = {rgb('FireBrick'), 1,1,1,1, rgb('DarkTurquoise')}

for i = [1 6]
opts.visCond = i; 
[popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
[eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
dataInPlots{1}=outDat{1};
eHandle{1}.CapSize =0;
title('min')
eHandle{1}.Color = clist{i};
end
% [num2str(mean(test1(opts.ensemblesToPlot)),2) ' ' num2str(mean(test2(opts.ensemblesToPlot)),2)]

%%
% %% Distance by Vis Response (will only work if consistent number of unique(visID)
% plotResponseByDistanceContrast(outVars,opts); %warning won't throw an error even if you have no contrasts
% %% Contrast Response Functions
% 
% opts.distBinCRF = [0 25 50:50:350]; %split Contrast Response Fun by distance
% opts.visAlphaCRF = 0.05; %visAlpha for looking just at vis responsive cells;
% 
% [outVars] = plotContrastResponseFunction(All,outVars,opts);
% 
% %% I think nothing below here...
% % %% Plot Mean Distance Responses
% 
% plotEnsembleDistanceResponse(outVars,100,1)
% 
% 
% 
% %% Plot Distance Plots
% 
% opts.distBins = 0:25:1000; %must be set to match popDist
% plotResponseByDistance(outVars,opts);
% 
% 
% %% Plot Distance with Filters
% 
% 
% opts.ensemblesToPlot = outVars.ensemblesToUse;
% opts.useVisCells =0;
% opts.useTunedCells =0;
% opts.minNumberOfCellsPerCondition = -1;
% sum(opts.ensemblesToPlot)
% opts.variableCellFun =  'outVars.pVisR{ind} > -0.05';
% [allPopVal] = subsetPopResponse(All,outVars,opts);
% 
% figure(105);clf;
% 
% ax1 =  subplot(1,2,1);
% s = scatter(outVars.ensMaxD(opts.ensemblesToPlot),allPopVal(opts.ensemblesToPlot),'filled');
% s.MarkerFaceColor = rgb('firebrick');
% hold on
% 
% xlabel(['Spread of Ensemble (\mum)'])
% ylabel('Population Mean Response')
% % title('OSIs by Ensemble Size')
% 
% r = refline(0);
% r.Color = rgb('grey');
% r.LineStyle =':';
% r.LineWidth = 2;
% 
% A = outVars.ensMaxD(opts.ensemblesToPlot)';
% B = allPopVal(opts.ensemblesToPlot)';
% [p Rsq, pVal] = simplifiedLinearRegression(A,B);
% 
% pl = plot(outVars.ensMaxD(opts.ensemblesToPlot)',p(1)*outVars.ensMaxD(opts.ensemblesToPlot)+p(2),'color',rgb('dimgrey'));
% title({'All Cells'; ['Fit pVal: ' num2str(pVal(1),3)]; ['R2: ' num2str(Rsq(1))] })
% 
% ax2 = subplot(1,2,2);
% 
% opts.variableCellFun =  'outVars.pVisR{ind} < 0.05';
% % opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & outVars.osi{ind}<0.25';
% 
% [visPopVal] = subsetPopResponse(All,outVars,opts);
% s = scatter(outVars.ensMaxD(opts.ensemblesToPlot),visPopVal(opts.ensemblesToPlot),'filled');
% s.MarkerFaceColor = rgb('firebrick');
% hold on
% 
% xlabel(['Spread of Ensemble (\mum)'])
% ylabel('Population Mean Response')
% % title('OSIs by Ensemble Size')
% 
% r = refline(0);
% r.Color = rgb('grey');
% r.LineStyle =':';
% r.LineWidth = 2;
% 
% A = outVars.ensMaxD(opts.ensemblesToPlot)';
% B = visPopVal(opts.ensemblesToPlot)';
% [p Rsq, pVal] = simplifiedLinearRegression(A,B);
% 
% pl = plot(outVars.ensMaxD(opts.ensemblesToPlot)',p(1)*outVars.ensMaxD(opts.ensemblesToPlot)+p(2),'color',rgb('dimgrey'));
% title({'Just Vis Cells'; ['Fit pVal: ' num2str(pVal(1),3)]; ['R2: ' num2str(Rsq(1))] })
% 
% linkaxes([ax1,ax2])
% 
% % %% Compare Distance responses
% % figure(102);clf
% % 
% % dataInPlots=[];
% % distTypes = {'min' 'geo' 'mean' 'harm' 'median' 'centroid'};
% % for i =[1 3]; %1:6
% %     disp(['working on ' distTypes{i}])
% %     opts.distType = distTypes{i}; %options: min geo mean harm
% %     opts.distBins = 0:10:350; %can be set variably 0:25:1000 is defaultt
% %     CellToUseVar = 'anal.cellsToInclude'; %[];
% %     [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
% %     ax = subplot(2,3,i);
% %     opts.distAxisRange = [0 350]; %[0 350] is stand
% %     [eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
% %     dataInPlots{i}=outDat{1};
% %     eHandle{1}.CapSize =0;
% %     title(distTypes{i})
% %     drawnow
% % end
% % disp('done')
% %%
% muPerPx = 800/512;
% opts.thisPlaneTolerance =15/muPerPx;% 15/muPerPx;
% opts.onePlaneTolerance = 30/muPerPx; %30/muPerPx;
% 
% recalcOffTargetRisk;
% %% Just a few with different binning
% figure(104);clf;
% dataInPlots =[];
% 
% ax = subplot(1,2,1);
% opts.distType = 'min';
% opts.distBins = 15:15:250; %can be set variably 0:25:1000 is defaultt
% opts.distAxisRange = [0 250]; %[0 350] is stand
% % CellToUseVar = 'anal.cellsToInclude & All(ind).out.anal.pVisR<0.05';%[];
% CellToUseVar = 'anal.cellsToInclude ';%[];
% 
% [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
% [eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
% dataInPlots{1}=outDat{1};
% eHandle{1}.CapSize =0;
% title('min')
% 
% ax = subplot(1,2,2);
% opts.distType = 'mean';
% opts.distBins = 15:15:500; %can be set variably 0:25:1000 is defaultt
% opts.distAxisRange = [0 450]; %[0 350] is stand
% CellToUseVar = 'anal.cellsToInclude';%[];
% [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
% [eHandle outDat] = plotDistRespGeneric(popRespDist,outVars,opts,ax);
% dataInPlots{2}=outDat{1};
% eHandle{1}.CapSize =0;
% title('mean')
% ylim([-0.075 0.075])
% disp('done')
% 
% %%
% %% Plot Split by Mean OSI
% opts.distType = 'min';
% opts.distBins =[0:25:150]; 
% opts.plotTraces = 0; 
% opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;% & outVars.ensOSI>0.3 & outVars.ensOSI<0.7;
% opts.criteriaToSplit = outVars.ensMaxD;
% opts.criteriaBins = [0 400 500 inf];
% % opts.useVisAndTunedCells =1; 
% opts.useVisCells = 1; 
% opts.useTunedCells = 1; 
% 
% plotDistByCriteria(All,outVars,opts,15)
% figure(16);
% ylim([-0.15 0.15]);
% 
% 
% %% Fail StimTest Vs Not Plots
% countUSC=[]
% for ind = 1:numExps
%     allStimmedCells = unique([All(ind).out.exp.holoTargets{:}]);
%     allStimmedCells(isnan(allStimmedCells)) = [];
%     tempCells = 1:size(All(ind).out.exp.zdfData,1);
%     asc = ismember(tempCells,allStimmedCells);
%     All(ind).out.anal.allStimmedCells = asc;
%     All(ind).out.anal.neverStimmedCells = ~asc;
% 
%      tc = ismember(tempCells,All(ind).out.exp.targetedCells);
%      All(ind).out.anal.allTargetedCells = tc;
%      All(ind).out.anal.neverTargetedCells = ~tc;
% 
%      if isfield(All(ind).out.exp,'holoRequest')
%          roiNaN =isnan(All(ind).out.exp.holoRequest.roiWeights);% | All(ind).out.exp.holoRequest.roiWeights>1.5;
%          unstimableCells = All(ind).out.exp.targetedCells(roiNaN);
%          usc = ismember(tempCells,unstimableCells);
%      else
%          usc = zeros(size(tempCells));
%      end
% 
%      countUSC(ind) = sum(usc);
%      disp(['Ind: ' num2str(ind) '. ' num2str(sum(usc)) ' cells unstimmable'])
% 
%      All(ind).out.anal.unstimableCells= usc;
% 
% end
% disp('made new stimmed vs neverstimmed cells')
% 
% %% Plot
% figure(104);clf
% ax = subplot(1,1,1);
% 
% opts.distBins = 0:5:350; %0:25:350; %can be set variably 0:25:1000 is defaultt
% opts.distType = 'min';
% opts.distAxisRange = [0 250]; %[0 350] is stand
% 
% backupEnsemblesToUse = outVars.ensemblesToUse;
% noUnstimableCount = find(countUSC==0);
% limEnsembleToUse = outVars.ensemblesToUse & ~ismember(outVars.ensIndNumber,noUnstimableCount);
% outVars.ensemblesToUse = limEnsembleToUse;
% disp(['Using only ' num2str(sum(outVars.ensemblesToUse)) ' Ensembles']);
% % 
% CellToUseVar = [];
% [popRespDistDefault] = popDistMaker(opts,All,CellToUseVar,0);
% p1 = plotDistRespGeneric(popRespDistDefault,outVars,opts,ax);
% p1{1}.Color=rgb('black');
% outVars.ensemblesToUse = backupEnsemblesToUse;
% hold on
% drawnow
% 
% % CellToUseVar = 'anal.neverTargetedCells';
% % [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
% % p2 = plotDistRespGeneric(popRespDist,outVars,opts,ax);
% % p2{1}.Color=rgb('Mediumblue');
% 
% drawnow
% 
% CellToUseVar = 'anal.unstimableCells';
% [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
% p3 = plotDistRespGeneric(popRespDist,outVars,opts,ax);
% p3{1}.Color=rgb('ForestGreen');
% % legend([p1{1} p2{1} p3{1}],'All Cells','Not Tested', 'Unstimable');
% legend([p1{1}  p3{1}],'All Cells', 'Unstimable');
% 
% ylim([-0.05 0.11])
% 
% tempDat1 = popRespDistDefault(limEnsembleToUse,1);
% tempDat2 = popRespDist(limEnsembleToUse,1);
% 
% ranksum(tempDat1,tempDat2)