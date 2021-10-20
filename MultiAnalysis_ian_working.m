%% Load Experiments/Setup
clear
close all


addpath(genpath('100spikesAnalysis'))
%% loadLists

oriLoadList;
% % allLoadList;

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

%exclude certain expression types:
uniqueExpressionTypes = outVars.uniqueExpressionTypes;
excludedTypes ={'AAV CamK2' 'Ai203' 'neo-IV Tre 2s' 'IUE CAG' };


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
    & ensembleOneSecond ... %cuts off a lot of the earlier
    & numCellsEachEns==10 ...
    & ensDate >= -210428 ...
    ...& outVars.hzEachEns == 10 ...
    ...& outVars.hzEachEns >= 9 & outVars.hzEachEns <= 12 ...
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

%% Create time series plot

[All, outVars] = createTSPlot(All,outVars);
% [All, outVars] = createTSPlotAllVis(All,outVars);

%% Optional group ensembles into small medium and large
numCellsEachEns = outVars.numCellsEachEnsBackup;
numCellsEachEns(numCellsEachEns < 10) = 3;
numCellsEachEns(numCellsEachEns > 12) = 20;

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
% [All, outVars] = createTSPlotByEnsSizeAllVis(All,outVars);

%% Distance Response Plots
opts.distBins = 0:25:1000; %must be set to match popDist
plotResponseByDistance(outVars,opts);


%% Compare Distance responses
figure(102);clf

distTypes = {'min' 'geo' 'mean' 'harm' 'median' 'centroid'};
for i =1:6
    disp(['working on ' distTypes{i}])
    opts.distType = distTypes{i}; %options: min geo mean harm
    opts.distBins = 0:1:350; %can be set variably 0:25:1000 is defaultt
    CellToUseVar = [];
    [popRespDist] = popDistMaker(opts,All,CellToUseVar,0);
    ax = subplot(2,3,i);
    opts.distAxisRange = [0 350]; %[0 350] is stand
    plotDistRespGeneric(popRespDist,outVars,opts,ax);
    title(distTypes{i})
    drawnow
end
disp('done')

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

opts.ensXaxis = 'osi'; % order, osi, dist, corr, size...
% plotCompareAllCellsTuning(outVars, opts);
opts.goodOSIthresh = 0.7; %ensemble OSI threshold
plotResponseByDifferenceinAnglePref(outVars, opts)

%% seperate by posNeg
[All outVars] = posNegIdentifiers(All,outVars,opts);
outVars = getRespByTuningDiffPosNeg(All, outVars);

plotResponseByDiffAnglePrefPosNeg(outVars,opts);

%% Pos vs Neg by Distance
opts.posNegThreshold = 0.0;
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
opts.distAxisRange = [0 350];

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

%%
opts.posNegThreshold = 0;%0.1;
[All outVars] = posNegIdentifiers(All,outVars,opts);


%% 2D Ensemble Response Plot
% opts.distType = 'min';
% [outVars] = grandDistanceMaker(opts,All,outVars);
% numEns = numel(outVars.ensStimScore);
% 
% distRange = [0 35];[50 150];%[50 150];
% for i = 1:numEns
%     responseToPlot(i) = nanmean(outVars.mRespEns{i}(...
%         outVars.distToEnsemble{i}>=distRange(1) & ...
%         outVars.distToEnsemble{i}<distRange(2) & ...
%         ~outVars.offTargetRiskEns{i} ...
%         ...& outVars.isRedByEns{i} ...
%         ));
% end
% % responseToPlot = outVars.popResponseEns; %outVars.mRespEns;
% yToPlot = outVars.ensOSI;
% xToPlot = outVars.ensMaxD;
% 
% ensemblesToUse = outVars.ensemblesToUse & outVars.numMatchedTargets>2 ;%& outVars.numMatchedTargets<12;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;% &  outVars.ensOSI>-0.25 & outVars.numMatchedTargets>=3 & outVars.ensMaxD<400;% outVars.numMatchedTargets>=3 &    ;
% 
% figure(142);clf
% scatter(xToPlot(ensemblesToUse),yToPlot(ensemblesToUse),100,responseToPlot(ensemblesToUse),'filled');
% colormap(puor)
% colorbar
% % caxis([-0.1 0.1])
% caxis([-0.25 0.25])
% xlabel({'Span of Ensemble (\mum)' ' (i.e. max distance between targets)'})
% ylabel('Ensemble OSI')


%% Plot Distance Curves in different Ori Bands

opts.distType = 'min';
opts.distBins =[15:15:150];% [0:25:400];
opts.distBins =[0:25:150];% [0:25:400];
% 
% opts.distType = 'harm'; 
% opts.distBins = [0:50:400];

[outVars] = grandDistanceMaker(opts,All,outVars);
numEns = numel(outVars.ensStimScore);


%this is where you change the criteria of what ensembles are included
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==1 &  outVars.ensOSI>0.7;% & outVars.meanEnsOSI>0.25;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse &  outVars.ensOSI>0.7 & outVars.meanEnsOSI>0.5;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI>0.7;% & outVars.meanEnsOSI>0.25;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;


oriVals = [NaN 0:45:315];
% numEns = numel(outVars.posCellbyInd);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
% 
    diffsPossible = [0 45 90 135 180];

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

       cellOris = oriVals(outVars.prefOris{ind});
       cellOrisDiff = abs(cellOris-outVars.ensPO(i));
       cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
%        
       cellOrisDiff(cellOrisDiff==135)=45;
       cellOrisDiff(cellOrisDiff==180)=0;


       %%This is where you change the criteria of what cells are included
    cellToUseVar = ~outVars.offTargetRiskEns{i}...
        & outVars.pVisR{ind} < 0.05 ...
        & outVars.osi{ind} > 0.25 ...
         ...& outVars.posCellbyInd{i} ... %if you want to only include cells that went up
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
opts.distAxisRange = [0 350];

figure(11);clf;hold on
figure(14);clf

% ax =subplot(1,1,1);
%  hold on
colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
dataForStats=[];
clear ax
for k = 1:numel(diffsPossible)
    figure(10);
    ax(k) =subplot(1,numel(diffsPossible),k);
    title(['Cells Pref Angle \Delta' num2str(diffsPossible(k)) '\circ'])
    [eHandle outData] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,ax(k));
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end
    eHandle{1}.CapSize =0;
    ylabel('Pop Response (Mean \DeltaF/F)')

            dataForStats(k,:) = outData{1}.dat(:,1);

    
    figure(11); tempax = subplot(1,1,1);
    [eHandle] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,tempax);
    eHandle{1}.Color = colorListOri{k};
    ylabel('Pop Response (Mean \DeltaF/F)')

    if k ==1 || k==3
        figure(14); tempax = subplot(1,1,1);
        [eHandle] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,tempax);
        eHandle{1}.Color = colorListOri{k};
        eHandle{1}.CapSize =0;
        ylabel('Pop Response (Mean \DeltaF/F)')

        if k==1
        delete(eHandle{end})
        end
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])

figure(10);
xlim([0 opts.distBins(end)])
ylim([-0.25 0.25])

figure(14)
xlim([0 opts.distBins(end)])
legend('Iso','Ortho')



figure(12);clf
datToPlot = squeeze(popToPlot(ensemblesToUse,1,:));
ciToPlot = nanstd(datToPlot,[],1)./sqrt(size(datToPlot,1))*1.93;
errorbar(nanmean(datToPlot),ciToPlot);
p = anova1(datToPlot,[],'off');
% p2 = signrank(datToPlot(:,1),datToPlot(:,3));
% p2 = signrank([datToPlot(:,1);datToPlot(:,5)],datToPlot(:,3));
% [~,p2] = ttest2(datToPlot(:,1),datToPlot(:,3));
% [~,p2] = ttest2([datToPlot(:,1);datToPlot(:,5)],datToPlot(:,3));


p2 = ranksum(squeeze(dataForStats(1,:)),squeeze(dataForStats(3,:)));
disp(['ranksum iso v ortho p = ' num2str(p2)]);

title({['Pvalue ' num2str(p)]...
    ['Iso vs Ortho PValue: ' num2str(p2)] })

% title('Cells by Tuning')
% ax2 =subplot(1,2,2);
% plotDistRespGeneric(popToPlotNeg,outVars,opts,ax2);
% title('Cells That go down')

%% Plot Distance Plots by Ensemble Tuning
opts.distType = 'min';
opts.distBins =[0:25:150]; [15:15:150];[10:20:150];% [0:25:400];

plotTraces=0;
%things to hold constant
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;% & outVars.meanEnsOSI>0.5 ;%  &  outVars.ensOSI>0.75;;
criteria = outVars.ensOSI;outVars.meanEnsOSI;
useableCriteria = criteria(ensemblesToUse);
bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins = [0 0.25 0.5 0.75 max(useableCriteria)];
bins = [0 0.33 0.7 max(useableCriteria)];
bins = [0 0.3 0.7 max(useableCriteria)];


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
            & outVars.osi{ind} > 0.25 ...
            ... & outVars.isRedByEns{i} ...
            ;
        
        popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')

figure(16);clf
figure(15);clf

colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
for k = 1:numel(bins)-1
    figure(16);ax(k) =subplot(1,numel(bins)-1,k);
    title(['Ens Criteria: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ])
    popToPlotTemp = popToPlot;
    popToPlotTemp(criteria<bins(k) | criteria>bins(k+1),:)=NaN;
    [eHandle outDat] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(k));
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end
    if plotTraces
        hold on
        distBins = opts.distBins;
        distBinSize = distBins(2)-distBins(1);
        plot(distBins(2:end)-distBinSize/2,outDat{1}.dat','color',rgb('grey'));
    end
    
    figure(15); tempax = subplot(1,1,1);
    [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
    eHandle{1}.Color = colorListOri{k};
end
linkaxes(ax)
figure(16);
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])
figure(15);
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])

%% Plot Distance Plots by Ensemble Tuning and Spacing
% opts.distType = 'harm'; 
% opts.distBins =[0:50:400];
opts.distType = 'min';
opts.distBins =[0:25:150];%[15:20:150];% [0:25:400];
opts.plotTraces =0;
% plotSpaceAndFeature(All,outVars,opts)

%things to hold constant
ensemblesToUse = ...
    outVars.ensemblesToUse ...
    & outVars.numCellsEachEnsBackup==10 ...
    ...& outVars.meanEnsOSI>0.5 ...
    ;
criteria =  outVars.ensMaxD; outVars.ensMeaD;%   outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse);
bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins = [0 400 500 inf];
% bins = [0 400  inf];
% bins = [0 400 500 inf];%bins = [0 200 250 inf];

% bins = [linspace(min(useableCriteria),max(useableCriteria),5)];
criteria2 =   outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria2(ensemblesToUse);
bins2 = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
% bins2 = [0 0.33 0.7 inf];
bins2 = [0 0.3 0.7 inf];
% bins2 = [0 0.25 0.5 inf];

% bins2 = [linspace(min(useableCriteria),max(useableCriteria),5)];

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

f = figure(555);
clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

    cellToUseVar = ~outVars.offTargetRiskEns{i}...
          & outVars.pVisR{ind} < 0.05 ...
          ... & outVars.osi{ind} > 0.25 ...
         ...& outVars.posCellbyInd{i} ...
        ;

        popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
close(f)
disp('Done')

figure(13);clf

colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
dataForStats =[];
c=0;
for i = 1:numel(bins2)-1
    for k = 1:numel(bins)-1
        c=c+1;
        
        ax(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
        popToPlotTemp = popToPlot;
        ensembleSelecter  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;
        popToPlotTemp(~ensembleSelecter,:)=NaN;
        [eHandle outDat] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(c));
        eHandle{1}.CapSize =0;
        if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
            eHandle{1}.Color = colorListOri{k};
        end
        dataForStats(i,k,:) = outDat{1}.dat(:,1);
        
        if opts.plotTraces
        hold on
        distBins = opts.distBins;
        distBinSize = distBins(2)-distBins(1);
        plot(distBins(2:end)-distBinSize/2,outDat{1}.dat','color',rgb('grey'));
        end
        
        title({...
            ['Ens Dist: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            ['Ens OSI: ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            ['Num Ens: ' num2str(sum(ensembleSelecter )) ] ...
            } )
        
        ylabel('Pop Response (Mean \DeltaF/F)')
        
        %     figure(11); tempax = subplot(1,1,1);
        %     [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
        %     eHandle{1}.Color = colorListOri{k};
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.25 0.3])

p1 = ranksum(squeeze(dataForStats(3,1,:)),squeeze(dataForStats(3,3,:)));
% [a p1] = ttest2(squeeze(dataForStats(3,1,:)),squeeze(dataForStats(3,3,:)));
disp(['Co Tuned Close vs Far p = ' num2str(p1)]);

p2 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(1,3,:)));
% [a p2] = ttest2(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(1,3,:)));
disp(['UnTuned Close vs Far p = ' num2str(p2)]);

p3 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(3,1,:)));
% [a p3] = ttest2(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(3,1,:)));
disp(['Close Tuned v Untuned p = ' num2str(p3)]);

p4 = ranksum(squeeze(dataForStats(1,3,:)),squeeze(dataForStats(3,3,:)));
% [a p4] = ttest2(squeeze(dataForStats(1,3,:)),squeeze(dataForStats(3,3,:)));
disp(['Far Tuned v Untuned p = ' num2str(p4)]);

% ylim([-1 1]);
% figure(11);
% xlim([0 opts.distBins(end)])
% ylim([-0.1 0.3])

%% Plot Ori Bands Close vs Far


opts.distType = 'min';
opts.distBins =[15:15:150];% [0:25:400];
opts.distBins =[0:25:150];% [0:25:400];
% 
% opts.distType = 'harm'; 
% opts.distBins = [0:50:400];

[outVars] = grandDistanceMaker(opts,All,outVars);
numEns = numel(outVars.ensStimScore);


%this is where you change the criteria of what ensembles are included
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 ...
    &  outVars.ensOSI>0.7;% & outVars.ensMaxD >400 & outVars.ensMaxD <500  ;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;


oriVals = [NaN 0:45:315];
% numEns = numel(outVars.posCellbyInd);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
% 
    diffsPossible = [0 45 90 135 180];

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

       cellOris = oriVals(outVars.prefOris{ind});
       cellOrisDiff = abs(cellOris-outVars.ensPO(i));
       cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
       
       cellOrisDiff(cellOrisDiff==135)=45;
       cellOrisDiff(cellOrisDiff==180)=0;


       %%This is where you change the criteria of what cells are included
    cellToUseVar = ~outVars.offTargetRiskEns{i}...
        & outVars.pVisR{ind} < 0.05 ...
        & outVars.osi{ind} > 0.25 ...
         ...& outVars.posCellbyInd{i} ... %if you want to only include cells that went up
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
opts.distAxisRange = [0 350];

figure(11);clf;hold on
figure(14);clf

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
    eHandle{1}.CapSize =0;
    ylabel('Pop Response (Mean \DeltaF/F)')

    figure(11); tempax = subplot(1,1,1);
    [eHandle] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,tempax);
    eHandle{1}.Color = colorListOri{k};
    ylabel('Pop Response (Mean \DeltaF/F)')

    if k ==1 || k==3
        figure(14); tempax = subplot(1,1,1);
        [eHandle outDat] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,tempax);
        eHandle{1}.Color = colorListOri{k};
        eHandle{1}.CapSize =0;
        ylabel('Pop Response (Mean \DeltaF/F)')

        if k==1
        delete(eHandle{end})
        end
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])

figure(10);
xlim([0 opts.distBins(end)])

figure(14)
xlim([0 opts.distBins(end)])
legend('Iso','Ortho')
ylim([-0.35 0.35])
axis square


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

%% Plot Distance Plots by Ensemble OSI and OSI of ensemble Members
%GREG LOOK AT THIS
% opts.distType = 'harm'; 
% opts.distBins =[0:50:400];
opts.distType = 'min';
opts.distBins =[0:25:150];[15:20:150];% [0:25:400];

plotTraces =0;
%things to hold constant
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 % & outVars.meanEnsOSI>0.5;
criteria =  outVars.meanEnsOSI;  outVars.ensMaxD; outVars.ensMeaD;%   outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse);
bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
% bins = [0  400 500 inf];
bins = [0 0.3 0.5 inf];
% bins = [0 400  inf];
% bins = [0 400 500 inf];%bins = [0 200 250 inf];

% bins = [linspace(min(useableCriteria),max(useableCriteria),5)];
criteria2 =   outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria2(ensemblesToUse);
bins2 = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
% bins2 = [0 0.33 0.7 inf];
bins2 = [0 0.3 0.7 inf];
%  bins2 = [0 0.25 0.5 inf];

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
dataForStats =[];
c=0;
for i = 1:numel(bins2)-1
    for k = 1:numel(bins)-1
        c=c+1;
        
        ax(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
        popToPlotTemp = popToPlot;
        ensembleSelecter  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;
        popToPlotTemp(~ensembleSelecter,:)=NaN;
        [eHandle outDat] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(c));
        eHandle{1}.CapSize =0;
        if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
            eHandle{1}.Color = colorListOri{k};
        end
        dataForStats(i,k,:) = outDat{1}.dat(:,1);
        
        if plotTraces
        hold on
        distBins = opts.distBins;
        distBinSize = distBins(2)-distBins(1);
        plot(distBins(2:end)-distBinSize/2,outDat{1}.dat','color',rgb('grey'));
        end
        
        title({...
            ['Mean OSI: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            ['Ens OSI: ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            ['Num Ens: ' num2str(sum(ensembleSelecter & ensemblesToUse)) ] ...
            } )
        
        ylabel('Pop Response (Mean \DeltaF/F)')
        
        %     figure(11); tempax = subplot(1,1,1);
        %     [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
        %     eHandle{1}.Color = colorListOri{k};
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.25 0.3])

p1 = ranksum(squeeze(dataForStats(3,1,:)),squeeze(dataForStats(3,3,:)));
disp(['Co Tuned Close vs Far p = ' num2str(p1)]);
p2 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(1,3,:)));
disp(['UnTuned Close vs Far p = ' num2str(p2)]);

p3 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(3,1,:)));
disp(['Close Tuned v Untuned p = ' num2str(p3)]);
p4 = ranksum(squeeze(dataForStats(1,3,:)),squeeze(dataForStats(3,3,:)));
disp(['Far Tuned v Untuned p = ' num2str(p4)]);
