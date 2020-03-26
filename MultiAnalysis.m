%% Load Experiments
addpath('SubFunctions')

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

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
opts.recWinRange = [0.5 1.5];% %from vis Start [1.25 2.5];


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

excludeInds = ismember(ensIndNumber,[3]); %Its possible that the visStimIDs got messed up




ensemblesToUse = numSpikesEachEns > 75 &...
    numSpikesEachEns <110 &...
    ...highVisPercentInd &...
    lowRunInds &...
    ensStimScore > 0.5 &... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    ~excludeInds &...
    numTrialsPerEns > 3;%10;%&...
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

%% Create time series plot

[All, outVars] = createTSPlot(All,outVars);

%% Optional group ensembles into small medium and large
numCellsEachEns = outVars.numCellsEachEnsBackup;
% numCellsEachEns(numCellsEachEns <=5) = 5;
% numCellsEachEns(numCellsEachEns > 10) = 20;

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
[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);




%% Red Cell Analysis

