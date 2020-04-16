%% Load Experiments/Setup
clear all
close all
% addpath('SubFunctions') %will, does this work for you?, nope....
% addpath('circStats')
% addpath('LoadLists')

%%

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

%%
% oriLoadList
% loadList = loadList(15);

% allLoadList;
SSTOriLoadList;
% loadPath = 'U:\ioldenburg\outputdata1'
% loadPath = 'C:\Users\ian\Dropbox\Adesnik\Data\outputdata1'
% loadPath = 'C:\Users\SabatiniLab\Dropbox\Adesnik\Data\outputdata1' %Ian Desktop
% loadPath = 'C:\Users\Will\Local Data\100spikes-results\outfiles-ori'
loadPath = 'E:\100spikes-results\outfiles-all'
%%
numExps = numel(loadList);
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
opts.stimsuccessZ = 0.3; %0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %fraction of ensemble that needs to be succsfull

%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded


[All, outVars] = cleanData(All,opts);

ensStimScore        = outVars.ensStimScore;
hzEachEns           = outVars.hzEachEns;
numCellsEachEns     = outVars.numCellsEachEns;
numSpikesEachStim   = outVars.numSpikesEachStim;
percentLowRunTrials = outVars.percentLowRunTrials;
numSpikesEachEns    = outVars.numSpikesEachEns;

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

[All, outVars] = meanMatrixVisandCorr(All,opts,outVars);

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

%% main Ensembles to Use section
% ensemblesToUse = numSpikesEachEns > 75 & numSpikesEachEns <125 & highVisPercentInd & ensIndNumber~=15 & ensIndNumber~=16; %& numCellsEachEns>10 ;

numTrialsPerEns =[];numTrialsPerEnsTotal=[]; 
for ind=1:numExps
    us=unique(All(ind).out.exp.stimID);
    
    for i=1:numel(us)
        trialsToUse = All(ind).out.exp.lowMotionTrials &...
            All(ind).out.exp.lowRunTrials &...
            All(ind).out.exp.stimSuccessTrial &...
            All(ind).out.exp.stimID == us(i) ;
        
        numTrialsPerEns(end+1)=sum(trialsToUse);
        numTrialsPerEnsTotal(end+1) = sum(All(ind).out.exp.stimID == us(i));
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
excludedTypes = {'AAV CamK2' 'Ai203'};
% excludedTypes = {'Ai203'};

exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);


%spot to add additional Exclusions
excludeInds = ismember(ensIndNumber,[]); %Its possible that the visStimIDs got messed up

%Options
opts.numSpikeToUseRange = [98 101];
opts.ensStimScoreThreshold = 0.5; % default 0.5
opts.numTrialsPerEnsThreshold = 3; 

ensemblesToUse = numSpikesEachEns > opts.numSpikeToUseRange(1) ...
    & numSpikesEachEns < opts.numSpikeToUseRange(2) ...
    ...& highVisPercentInd ...
    & lowRunInds ...
    & ensStimScore > opts.ensStimScoreThreshold ... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    & ~excludeInds ...
    & numTrialsPerEns > opts.numTrialsPerEnsThreshold ... ;%10;%&...
    & ~ensHasRed ...
    & ~excludeExpressionType ...
    & ~ensMissedTarget ...
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

disp(['Fraction of Ens correct Size: ' num2str(mean(numSpikesEachEns > 75 & numSpikesEachEns <110))]);
disp(['Fraction of Ens highVis: ' num2str(mean(highVisPercentInd))]);
disp(['Fraction of Ens lowRun: ' num2str(mean(lowRunInds))]);
disp(['Fraction of Ens high stimScore: ' num2str(mean(ensStimScore>opts.ensStimScoreThreshold))]);
disp(['Fraction of Ens high trial count: ' num2str(mean(numTrialsPerEns>opts.numTrialsPerEnsThreshold))]);
disp(['Fraction of Ens No ''red'' cells shot: ' num2str(mean(~ensHasRed))]);
disp(['Fraction of Ens usable Expression Type: ' num2str(mean(~excludeExpressionType))]);
disp(['Fraction of Ens enough targets detected by s2p: ' num2str(mean(~ensMissedTarget))]);

disp(['Total Fraction of Ens Used: ' num2str(mean(ensemblesToUse))]);

%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end
%% Create time series plot

[All, outVars] = createTSPlot(All,outVars);

%% Optional group ensembles into small medium and large
numCellsEachEns = outVars.numCellsEachEnsBackup;
numCellsEachEns(numCellsEachEns <=5) = 5;
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
plotResponseBySize(outVars)
plotPopResponseBySession(All,outVars)
plotPopResponseByExpressionType(All,outVars);
[All outVars] = createTSPlotByEnsSize(All,outVars);
%% Distance Response Plots
plotResponseByDistance(outVars,opts);
plotResponseByDistanceContrast(outVars,opts); %warning won't throw an error even if you have no contrasts
%% Contrast Response Functions

opts.distBinCRF = [50:25:350]; %split Contrast Response Fun by distance
opts.visAlphaCRF = 10.05; %visAlpha for looking just at vis responsive cells;

[outVars] = plotContrastResponseFunction(All,outVars,opts);

%% Orientation Tuning and OSI
% leaving this a mess here for now to compare the different prefOri
% calculations.... will tidy up later
[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)

% compare max vs circular mean for determining PO
po = outVars.prefOris{1};
oris = 0:45:315;

osi = outVars.osi{1};
osi(po==1)=[];
po(po==1)=[];

podeg = oris(po-1);

figure(1)
clf
subplot(1,3,1)
scatter(podeg, outVars.circTuning{1}, [], osi, 'filled')
ylabel('Circular PO')
xlabel('Max PO')
c = colorbar;
c.Label.String = 'OSI';

subplot(1,3,2)
scatter(podeg, outVars.circTuning{1}, [], outVars.circVar{1}, 'filled')
ylabel('Circular PO')
xlabel('Max PO')
c = colorbar;
c.Label.String = 'Circ. Var (deg)';

% do OSI on new fancy tuning method
prefOri = outVars.circTuning{1};
ortho1 = mod(prefOri - 90, 135);
ortho2 = mod(prefOri + 90, 135);
orthoOri = cat(1,ortho1, ortho2);
% but have to cast back to 0:45:135
prefOriBinned = interp1(oris, oris, prefOri, 'nearest', 'extrap');
orthoOriBinned = interp1(oris, oris, orthoOri, 'nearest', 'extrap');
% then go back to visID to get idxs
for i=1:numel(prefOriBinned)
    o = prefOriBinned(i);
    oo = orthoOriBinned(:,i);
    prefIDs(i) = find(o==oris);
    orthoIDs1(i) = find(oo(1,:)==oris);
    orthoIDs2(i) = find(oo(2,:)==oris);
end
orthoIDs = cat(1,orthoIDs1, orthoIDs2);

%%
oriCurveBL = outVars.circCurves{1}';
% oriCurveBL = curves - min(curves);
OSI=[];
for i=1:numel(prefOriBinned)
    OSI(i) = (oriCurveBL(prefIDs(i),i) - mean(oriCurveBL(orthoIDs(:,i)',i)))...
        / (oriCurveBL(prefIDs(i),i) + mean(oriCurveBL(orthoIDs(:,i)',i)));
end

subplot(1,3,3)
scatter(podeg, outVars.circTuning{1}, [], OSI, 'filled')
ylabel('Circular PO')
xlabel('Max PO')
c = colorbar;
c.Label.String = 'OSI by circ tuning';
%% Red Cell Analysis (will only run if you have the red section on all your recordings). 
[outVars] = plotResponseOfRedCells(All,outVars,opts);


%% Red Distance section
for ind = 1:numExps;
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

%% Pos vs Neg 
[All outVars] = posNegIdentifiers(All,outVars,opts);
opts.distType = 'min';

numEns = numel(outVars.posCellbyInd);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
    cellToUseVar = outVars.posCellbyInd{i};
    popToPlotPos(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

     cellToUseVar = outVars.negCellbyInd{i};
    popToPlotNeg(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

end
disp('Done')

%%
figure(10);clf
ax =subplot(1,2,1);
plotDistRespGeneric(popToPlotPos,outVars,opts,ax);
title('Cells that go up')
ax2 =subplot(1,2,2);
plotDistRespGeneric(popToPlotNeg,outVars,opts,ax2);
title('Cells That go down')

%% Correlation Pick One. Option A. Vis activity from interleaved Trials
%Functions are Mutually Exclusive.
[All, outVars] = defineCorrelationTypes(All,outVars); %Caution this and below are mutually exclusive
plotEnsembleCorrelationResponse(outVars,200,1);

opts.CorrToPlot = 'AllCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'
[outVars] = plotCorrelationResponse(All,outVars,opts.CorrToPlot);
%% Correlation Pick One. Option B. Vis Activity from out.vis epoch
[All, outVars] = defineCorrelationTypesOnVis(All, outVars); %Caution this and above are mutually exclusive
plotEnsembleCorrelationResponse(outVars,300,1)

opts.CorrToPlot = 'AllCorr'; % Options are: 'SpontCorr' 'AllCorr' AllMCorr' 'SignalCorr' and 'NoiseCorr'
[outVars] = plotCorrelationResponse(All,outVars,opts.CorrToPlot);

%% Distance same idea as above
[All, outVars] = defineDistanceTypes(All, outVars);
plotEnsembleDistanceResponse(outVars,100,1)

%%

