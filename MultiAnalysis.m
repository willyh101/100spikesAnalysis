%% Load Experiments/Setup
clear all
close all
addpath('SubFunctions')
addpath('circStats')


[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

%%
oriLoadList
loadList = loadList(15);
% loadPath = 'U:\ioldenburg\outputdata1'
% loadPath = 'C:\Users\ian\Dropbox\Adesnik\Data\outputdata1'
% loadPath = 'C:\Users\SabatiniLab\Dropbox\Adesnik\Data\outputdata1' %Ian Desktop
% loadPath = 'C:\Users\Will\Local Data\100spikes-results\outfiles-ori'
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
%% clean Data

opts.FRDefault=6;
opts.recWinRange = [0.5 1.5];% %from vis Start in s [1.25 2.5];


%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %over this number is a succesfull stim
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



%% Make all dataPlots into matrixes of mean responses
%%Determine Vis Responsive and Process Correlation

opts.visAlpha = 0.05;

%oftarget risk params
opts.thisPlaneTolerance = 11.25;%7.5;%1FWHM%10; %in um;% pixels
opts.onePlaneTolerance = 22.5;%15;%2FWHM %20;
opts.distBins =  [0:25:1000];

[All, outVars] = meanMatrixVisandCorr(All,opts,outVars);

visPercent = outVars.visPercent;
ensIndNumber =outVars.ensIndNumber;


%% Optional: Calc pVisR from Visual Epoch [CAUTION: OVERWRITES PREVIOUS pVisR]
[All, outVars] = CalcPVisRFromVis(All,opts,outVars);
visPercent = outVars.visPercent;

%% if there is a red section
[outVars] = detectShotRedCells(All,outVars)
ensHasRed = outVars.ensHasRed;
%% main Ensembles to Use section
% ensemblesToUse = numSpikesEachEns > 75 & numSpikesEachEns <125 & highVisPercentInd & ensIndNumber~=15 & ensIndNumber~=16; %& numCellsEachEns>10 ;

numTrialsPerEns =[];
for ind=1:numExps
    us=unique(All(ind).out.exp.stimID);
    
    for i=1:numel(us)
        trialsToUse = All(ind).out.exp.lowMotionTrials &...
            All(ind).out.exp.lowRunTrials &...
            All(ind).out.exp.stimSuccessTrial &...
            All(ind).out.exp.stimID == us(i) ;
        
        numTrialsPerEns(end+1)=sum(trialsToUse);
    end
    
end
numTrialsPerEns(numSpikesEachStim==0)=[];


highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<0.05)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));

excludeInds = ismember(ensIndNumber,[]); %Its possible that the visStimIDs got messed up




ensemblesToUse = numSpikesEachEns > 75 ...
    & numSpikesEachEns <110 ...
    ... & highVisPercentInd ...
    & lowRunInds ...
    & ensStimScore > 0.5 ... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    & ~excludeInds ...
    & numTrialsPerEns > 10 ... ;%10;%&...
    & ~ensHasRed ...
    ;
%& numCellsEachEns>10 ;

indsSub = ensIndNumber(ensemblesToUse);
IndsUsed = unique(ensIndNumber(ensemblesToUse));

sum(ensemblesToUse)

outVars.ensemblesToUse      = ensemblesToUse;
outVars.IndsUsed            = IndsUsed;
outVars.indsSub             = indsSub;
outVars.numTrialsPerEns     = numTrialsPerEns;
outVar.highVisPercentInd    = highVisPercentInd;
outVar.lowRunInds           = lowRunInds;

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

%% Optional Reset ensemble grouping to default
numCellsEachEns = outVars.numCellsEachEnsBackup;
outVars.numCellsEachEns= numCellsEachEns;

%% Basic Response Plots
plotAllEnsResponse(outVars)
plotResponseBySize(outVars)
plotPopResponseBySession(All,outVars)
%% Distance Response Plots
plotResponseByDistance(outVars,opts);
plotResponseByDistanceContrast(outVars,opts); %warning won't throw an error even if you have no contrasts
%% Contrast Response Functions

%eh I'll do it later

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
%% Red Cell Analysis
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

