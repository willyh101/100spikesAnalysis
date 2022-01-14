% Greg data selector ver2
% this will allow you to select data from the experiment, looking at
% responses to ensembles based off of tuning of the stimmed ensemble,
% currently written to look at one experiment at a time

%% setup junk
clear
close all

% addpath(genpath('100spikesAnalysis'))

restoredefaultpath;
folder = fileparts(which('GregAnalysis_v4.m')); 
addpath(genpath(folder));
rmpath(folder)

%%

load('/Users/gregoryhandy/Research_Local/outputdata1/210903_I151_3_outfile.mat')

%% load data

% basePath = 'T:\Outfiles'; % where ever your files are
basePath = '/Users/gregoryhandy/Research_Local/outputdata1';

[loadList, loadPath ]= uigetfile(basePath,'MultiSelect','on');
addpath(genpath(loadPath))

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

%%

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
excludedTypes ={'AAV CamK2' 'Ai203' 'neo-IV Tre 2s' 'IUE CAG' };%'SepW1 CAG 2s'};


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
    ...& ensDate >= -210428 ...
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

%% Orientation Tuning and OSI

[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
[All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically 
[All, outVars] = defineDistanceTypes(All, outVars);

All(ind).out.exp.ensMaxD =  outVars.ensMaxD;

%% the new thing for greg
% this will get you responses from cells in the experiment epoch

% for now we hard-code what ind we are looking at
% since there is just 1, set it to 1
% in the future, this could be all put into a for loop
ind = 1;

% select your ensemble tuning
% to view possible tunings to use run this: 
% All(ind).out.exp.ensPO

% choose the ensemble number you want (index into All(ind).out.exp.ensPO)

%  [NaN     0    45    90   135   180   225   270   315]
oriToUse = [45 135];

ensIDs = [find(All(ind).out.exp.ensPO==oriToUse(1)) find(All(ind).out.exp.ensPO==oriToUse(2))];
ensOri = [oriToUse(1)*ones(size(find(All(ind).out.exp.ensPO==oriToUse(1)),2),1);...
    oriToUse(2)*ones(size(find(All(ind).out.exp.ensPO==oriToUse(2)),2),1)];

ensOri'
All(ind).out.exp.ensMaxD(ensIDs)

ensCoTuned = All(ind).out.exp.meanEnsOSI(ensIDs)>0.2 & All(ind).out.exp.ensOSI(ensIDs)>0.2;

% 
% ensIDs = ensIDs(8);
% ensOri = ensOri(8);

% Below chance: ensOri(4) at 10%
% Above chance: ensIDs(5), ensIDs(6), ensIDs(7), ensIDs(8)

% ensIDs = ensIDs(All(ind).out.exp.ensMaxD(ensIDs)>590);
% ensOri = ensOri(All(ind).out.exp.ensMaxD(ensIDs)>590);

%%
% ensIDs = [ensIDs(5) ensIDs(7)];
% ensOri = [ensOri(5) ensOri(7)];
%%
trialsToUse_matrix=[];
orisTrialsMatrix = [];
cellsToUse_matrix = [];
ensMaxD_trials = [];
ensCoTuned_trials = [];
for ii = 1:length(ensIDs)
    ens = ensIDs(ii);
    
    s = outVars.ensHNumber(ens);
    us = unique(All(ind).out.exp.stimID);
    u = us(s);
    
    % trials to use
    trialsToUse_matrix(ii,:) = ismember(All(ind).out.exp.stimID, u) &...
        All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        ismember(All(ind).out.exp.visID, [0 1]);
    
    ensMaxD_trials(ii,:) = trialsToUse_matrix(ii,:)*All(ind).out.exp.ensMaxD(ensIDs(ii));
    ensCoTuned_trials(ii,:) = trialsToUse_matrix(ii,:)*ensCoTuned(ii);
    
    orisTrialsMatrix(ii,:) = ensOri(ii)*trialsToUse_matrix(ii,:);
    % cells to use
    holo = All(ind).out.exp.stimParams.roi{s};  % tracks stimID -> hologram
    tg = All(ind).out.exp.rois{holo}; % this is cells in ensemble
    offTargetRisks = All(ind).out.anal.offTargetRisk(holo,:);
    ROIinArtifact  = All(ind).out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
    visCells = All(ind).out.anal.pVisR < 0.05; % visually responsive
    orisUsed = [nan 0:45:315];
    prefOris = arrayfun(@(x) orisUsed(x), All(ind).out.anal.prefOri);
    
    % comment out parts if you want to change subsets of cells
    % note ~offTargetRisks excludes the targeted cells, if you want to keep the
    % off targets and exclude the targets (for this specific hologram) you can
    % use: ~ismember(1:numel(offTargetRisks), tg)
    cellsToUse_matrix(ii,:) = ~offTargetRisks &...
        ~ROIinArtifact' &...
        visCells;
%     cellsToUse_matrix(ii,:) = ~offTargetRisks;

    
    disp('....')
    disp(['Ensemble is tuned for ' num2str(All(ind).out.exp.ensPO(ens))])
    disp(['Ensemble OSI: ' num2str(All(ind).out.exp.ensOSI(ens))])
    disp(['Mean OSI: ' num2str(All(ind).out.exp.meanEnsOSI(ens))])
    disp(['Ensemble has a stim score of ' num2str(All(ind).out.exp.ensStimScore(ens))])
end

if length(ensIDs) > 1
    cellsToUse = logical(prod(cellsToUse_matrix));
    holoTrialsToUse = logical(sum(trialsToUse_matrix));
    holoOrisTrials = sum(orisTrialsMatrix);
    holoOrisTrials = holoOrisTrials(holoTrialsToUse);
else
    cellsToUse = logical(cellsToUse_matrix);
    holoTrialsToUse = logical(trialsToUse_matrix);
    holoOrisTrials = orisTrialsMatrix;
    holoOrisTrials = holoOrisTrials(holoTrialsToUse);
end

disp(['Total of ' num2str(sum(holoTrialsToUse)) ' trials used.'])
disp(['Total of ' num2str(sum(cellsToUse)) ' cells used.'])

% here are you traces from the experiment
% cells x time x trials
% note: these are not trialwise baselined (but could be)
traces = All(ind).out.exp.dataToUse(cellsToUse, :, holoTrialsToUse);




ensMaxD_trials = sum(ensMaxD_trials(:,holoTrialsToUse));
ensCoTuned_trials = sum(ensCoTuned_trials(:,holoTrialsToUse));

%% vis epoch
% if you want the responses from the vis epoch, for the same cells, run
% this code

% you can change to out.vis.dfData if it exists
visData = All(ind).out.vis.zdfData;

% trials to use
vs = All(ind).out.vis.visID;
orisUsed = [nan 0:45:315];
oris = arrayfun(@(x) orisUsed(x), vs);

% select which oris to look at here
orisToUseVis = oriToUse; % copy from exp  section

trialsToUse = ismember(oris, orisToUseVis);

% cells to use is copied from above
% result is cells x time x trials
% note: these are now trialwise baselined (but could be)
tracesVis = visData(cellsToUse, :, trialsToUse);
trialAngle = oris(trialsToUse);

%%

num_training_trials = min(floor(sum(ismember(oris, orisToUseVis(1)))/2),...
    floor(sum(ismember(oris, orisToUseVis(2)))/2));
num_testing_trials = num_training_trials;
num_neurons = size(tracesVis,1);

%%

opts.recWinRange = [0.5 1.5];
% timeframes for the visual stimulus 
recWinSec=opts.recWinRange+All(ind).out.vis.visStart;
winToUse = round(recWinSec*All(ind).out.info.FR);
% winToUse = [1 size(tracesVis,2)];
bwinToUse = max(round([0 All(ind).out.vis.visStart]*All(ind).out.info.FR),[1 1]);

%%

vis_ave = squeeze(mean(tracesVis(:,winToUse(1):winToUse(2),:),2));
training_matrix = [];
trainingLabels = [];

testing_matrix=[];
testingLabels=[];
for ii = 1:length(orisToUseVis)
    trainingTrials = find(trialAngle==orisToUseVis(ii),num_training_trials);
    training_matrix = [training_matrix; vis_ave(:,trainingTrials)'];
    trainingLabels = [trainingLabels; orisToUseVis(ii)*ones(num_training_trials,1)];
    
    testingTrials = find(trialAngle==orisToUseVis(ii),num_testing_trials,'last');
    testing_matrix = [testing_matrix; vis_ave(:,testingTrials)'];
    testingLabels = [testingLabels; orisToUseVis(ii)*ones(num_testing_trials,1)];
end

%% Train 200 classifiers to assess success of classifying holo stim tuning
num_trainings = 200;
overall = zeros(1,num_trainings);
pref = zeros(1,num_trainings);
ortho = zeros(1,num_trainings);
trial_correct = [];
for training_loop = 1:num_trainings

    %% Fit the logistic classifier on the training set of visual stims
    log_classifier =fitclinear(training_matrix,trainingLabels,'Learner','logistic','Regularization','lasso');
    
    %% Using the classifier, predict the testing set
    predictions = log_classifier.predict(testing_matrix);
    correct = 0;
    for ii = 1:length(predictions)
        if predictions(ii) == orisToUseVis(1) && ii <= num_testing_trials
            correct = correct + 1;
        elseif predictions(ii) == orisToUseVis(2) && ii > num_testing_trials
            correct = correct + 1;
        end
    end
    
    % It's usually 100% accurate. If it isn't, print out the accuracy
    if round(correct/length(predictions)*100) ~=100
        fprintf('Testing error in visual stim is %.2f',round(correct/length(predictions)*100));
    end
    
    %% Use the classifier to predit holographic tuning 
    holo_ave = squeeze(mean(traces(:,winToUse(1):winToUse(2),:),2))';
    holo_Predictions = log_classifier.predict(holo_ave);
    
    holoCorrect = 0;
    prefCorrect = 0;
    orthoCorrect = 0;
    for ii = 1:length(holo_Predictions)
        if holo_Predictions(ii) == holoOrisTrials(ii)
            trial_correct(ii,training_loop) = 1;
            holoCorrect = holoCorrect + 1;
            if holo_Predictions(ii) == oriToUse(1)
                prefCorrect =  prefCorrect + 1;
            else
                orthoCorrect =  orthoCorrect + 1;
            end
        end
    end
    overall(training_loop) = holoCorrect/length(holo_Predictions)*100;
    pref(training_loop) = prefCorrect/length(find(holoOrisTrials==oriToUse(1)))*100;
    ortho(training_loop) = orthoCorrect/length(find(holoOrisTrials==oriToUse(2)))*100;

% fprintf('Classifier overall percent correct: %.2f \n',overall)
% fprintf('Classifier pref percent correct: %.2f \n',pref)
% fprintf('Classifier ortho percent correct: %.2f \n',ortho)
% 
% length(find(holoOrisTrials==180))/length(holoOrisTrials)

end

%%

figure(32)
boxplot([overall,pref,ortho],...
    [ones(1,length(overall)),2*ones(1,length(overall)),3*ones(1,length(overall))])

%%
class_trial_correct = [];
for ii =1:length(ensIDs)
temp = (ensMaxD_trials == All(ind).out.exp.ensMaxD(ensIDs(ii)));
class_trial_correct(ii) = mean(sum(trial_correct(temp,:),2)/num_trainings);
end
plot(All(ind).out.exp.ensMaxD(ensIDs),class_trial_correct,'.','markersize',16)

