%%
% Trains (and tests) a binary classifier on neuronal activity evoked by a
% visual stimulus, and then attempts to classify the ensemble orientation preference
% based on the activity evoked from holographic stimulation
%
% Options
%   classLearner: '', 'svm', 'logistic'
%   classReg: '', 'lasso', 'ridge'
%%
function [holoTestAcc,visTestAcc,holoEnsTestAcc,holoEnsCotuned,ensIDs]=...
    binary_classifier_vis(outVars,All,oriToUse,classLearner,classReg,includedCells)

% ensIDs = [find(All.out.exp.ensPO==oriToUse(1) & outVars.ensemblesToUse==1)...
%     find(All.out.exp.ensPO==oriToUse(2) & outVars.ensemblesToUse==1)];
% ensOri = [oriToUse(1)*ones(size(find(All.out.exp.ensPO(outVars.ensemblesToUse)==oriToUse(1)),2),1);...
%     oriToUse(2)*ones(size(find(All.out.exp.ensPO(outVars.ensemblesToUse)==oriToUse(2)),2),1)];


ensIDs = find((All.out.exp.ensPO==oriToUse(1) | All.out.exp.ensPO==oriToUse(2))  ...
    & outVars.ensemblesToUse==1);



%% Preallocation

numTotalTrials = length(All.out.exp.lowMotionTrials);
numTotalCells = length(All.out.exp.allDepth);

trialsToUse_matrix=zeros(length(ensIDs), numTotalTrials);
orisTrialsMatrix = zeros(length(ensIDs), numTotalTrials);
ensMaxD_trials = zeros(length(ensIDs), numTotalTrials);
cellsToUse_matrix = zeros(length(ensIDs), numTotalCells);
tg_matrix = zeros(length(ensIDs), numTotalCells);

%% Loop through and store the ensemble data
for ii = 1:length(ensIDs)
    
    ensID = ensIDs(ii);
    
    s = outVars.ensHNumber(ensID);
    us = unique(All.out.exp.stimID);
    u = us(s);
    
    % Trials to use
    trialsToUse_matrix(ii,:) = ismember(All.out.exp.stimID, u) &...
        All.out.exp.lowMotionTrials &...
        All.out.exp.lowRunTrials &...
        All.out.exp.stimSuccessTrial &...
        ismember(All.out.exp.visID, [0 1]);
    
    % Trial information
    ensMaxD_trials(ii,:) = trialsToUse_matrix(ii,:)*All.out.exp.ensMaxD(ensID);
    orisTrialsMatrix(ii,:) = trialsToUse_matrix(ii,:)*All.out.exp.ensPO(ensID);
    
    %% Cells to use (depends on 'includedCells' input variable)
    holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
    % tg = All.out.exp.rois{holo}; % this is cells in ensemble
    % fixed finding the cells in the ensembles?
    tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));
    
    offTargetRisks = All.out.anal.offTargetRisk(holo,:);
    ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
    visCells = All.out.anal.pVisR < 0.05; % visually responsive
    
    % The same ensemble can be used under different stim conditions,
    % but that is accounted for above, so that the right trials are
    % used. Here is just a quick check that this is no mismatch in
    % targets
    % 22/02/09: ens~=holo only occurs in 211021_W40_2_outfile.mat
    if ensID~=holo && sum(All.out.exp.holoTargets{ensID}-All.out.exp.holoTargets{holo})>0
        error('Something went wrong the the ensemble/holo number')
    end
    
    % Keep the wanted cells from this ensemble
    % Note (Will): ~offTargetRisks excludes the targeted cells (by default),
    % If you want to keep the off targets and exclude the targets
    % you can use ~ismember(1:numel(offTargetRisks), tg)
    if strcmp(includedCells,'default')
        cellsToUse_matrix(ii,:) = ~offTargetRisks &...
            ~ROIinArtifact' &...
            visCells;
    elseif strcmp(includedCells,'all')
        cellsToUse_matrix(ii,:) = ones(size(offTargetRisks));
    elseif strcmp(includedCells,'notStimmedOrOT')
        cellsToUse_matrix(ii,:) = ~offTargetRisks;
    elseif strcmp(includedCells,'notStimmed')
        cellsToUse_matrix(ii,:) = ~ismember(1:numel(offTargetRisks), tg_adj);
    elseif strcmp(includedCells,'OT')
        tempMatrix = zeros(size(offTargetRisks));
        tempMatrix(tg_adj) = 1;
        tg_matrix(ii,:) = tempMatrix;
        cellsToUse_matrix(ii,:) = offTargetRisks;
    elseif strcmp(includedCells,'stimmedAndOT')
        cellsToUse_matrix(ii,:) = offTargetRisks;
    elseif strcmp(includedCells,'stimmed')
        cellsToUse_matrix(ii,:) = ismember(1:numel(offTargetRisks), tg_adj);
    else
        error('Unknown cell selection');
    end
end


%% Keep only the needed cells 
if strcmp(includedCells,'stimmedAndOT') || strcmp(includedCells,'stimmed')
    cellsToUse = sum(cellsToUse_matrix,1)>0;
elseif strcmp(includedCells,'OT') %special case of excluding stimmed cell
    tg_matrix = sum(tg_matrix,1)>0;
    cellsToUse_matrix = sum(cellsToUse_matrix,1)>0;
    cellsToUse = (cellsToUse_matrix-tg_matrix)>0;
else % only include cells that were never stimmed/off-targets
    cellsToUse = logical(prod(cellsToUse_matrix));
end
holoTrialsToUse = logical(sum(trialsToUse_matrix,1));
holoOrisTrials = sum(orisTrialsMatrix,1);
holoOrisTrials = holoOrisTrials(holoTrialsToUse);
ensMaxD_trials = sum(ensMaxD_trials(:,holoTrialsToUse),1);

% Store the zdfData from the right cells and trials 
tracesHolo = All.out.exp.zdfData(cellsToUse, :, holoTrialsToUse);

%% Now get the data from vis epochs

visData = All.out.vis.zdfData;

% Figure out which vis trials to use
vs = All.out.vis.visID;
vs(vs==0) = []; % Band-aid on an unknown issue; vs should be between 1 and 9
orisUsed = [nan 0:45:315];
oris = arrayfun(@(x) orisUsed(x), vs);
% select which oris to look at here
orisToUseVis = oriToUse; 
trialsToUse = ismember(oris, orisToUseVis);

% cells to use is copied from above (e.g., ~offTargets)
% result is cells x time x trials
tracesVis = visData(cellsToUse, :, trialsToUse);
trialAngle = oris(trialsToUse);

%% Optional pre-processing: zero-mean before stimulus across trials
% baselineActVis = squeeze(mean(mean(tracesVis(:,1:6,:),2)));
% for ii = 1: size(tracesVis,3)
%     tracesVis(:,:,ii) = tracesVis(:,:,ii)-baselineActVis(ii);
% end
% 
% baselineActHolo = squeeze(mean(mean(tracesHolo(:,1:6,:),2)));
% for ii = 1: size(tracesHolo,3)
%     tracesHolo(:,:,ii) = tracesHolo(:,:,ii)-baselineActHolo(ii);
% end

%% Create the training and testing matrices
num_training_trials = min(floor(sum(ismember(oris, orisToUseVis(1)))/2),...
    floor(sum(ismember(oris, orisToUseVis(2)))/2));
num_testing_trials = num_training_trials;

% timeframes for the visual stimulus
opts.recWinRange = [0.5 1.5];
recWinSec=opts.recWinRange+All.out.vis.visStart;
winToUse = round(recWinSec*All.out.info.FR); % Use the activity when the visual stim was present
% winToUse = [1 size(tracesVis,2)]; % Alternative; use all data

visAve = squeeze(mean(tracesVis(:,winToUse(1):winToUse(2),:),2));
size(visAve)
trainingMatrix = []; trainingLabels = [];
testingMatrix=[]; testingLabels=[];
for ii = 1:length(orisToUseVis)
    trainingTrials = find(trialAngle==orisToUseVis(ii),num_training_trials);
    trainingMatrix = [trainingMatrix; visAve(:,trainingTrials)'];
    trainingLabels = [trainingLabels; orisToUseVis(ii)*ones(num_training_trials,1)];
    
    testingTrials = find(trialAngle==orisToUseVis(ii),num_testing_trials,'last');
    testingMatrix = [testingMatrix; visAve(:,testingTrials)'];
    testingLabels = [testingLabels; orisToUseVis(ii)*ones(num_testing_trials,1)];
end

%% Train the binary classifier on the data from the visual stims
binaryClassifier =fitclinear(trainingMatrix,trainingLabels,'Learner',...
    classLearner,'Regularization',classReg);

%% Predict the testing set with classifier
visPredictions = binaryClassifier.predict(testingMatrix);
visCorrect = 0;
for ii = 1:length(visPredictions)
    if visPredictions(ii) == orisToUseVis(1) && ii <= num_testing_trials
        visCorrect = visCorrect + 1;
    elseif visPredictions(ii) == orisToUseVis(2) && ii > num_testing_trials
        visCorrect = visCorrect + 1;
    end
end

% Save the testing accuracy for the classifier
visTestAcc = round(visCorrect/length(visPredictions)*100);

%% Predict holographic tuning with classifier
holoAve = squeeze(mean(tracesHolo(:,winToUse(1):winToUse(2),:),2))';
holoPredictions = binaryClassifier.predict(holoAve);

%% Store the overall performance across all ensembles
holoCorrect = (holoPredictions' ==holoOrisTrials);
holoTestAcc = sum(holoCorrect)/length(holoCorrect)*100;

%% Store performance for individual ensembles
holoEnsTestAcc = [];
holoEnsCotuned = [];
% Only save this data if the classifier was accurate on vis stimuli
if visTestAcc>=70
    for ii =1:length(ensIDs)
        
        if All.out.exp.meanEnsOSI(ensIDs(ii))>0.5 && All.out.exp.ensOSI(ensIDs(ii))>0.5
            holoEnsCotuned(ii) = 1;
        else
            holoEnsCotuned(ii) = 0;
        end
        
        % 'clever' way to find the ensIDs in the trial list
        temp = (ensMaxD_trials == All.out.exp.ensMaxD(ensIDs(ii)));
        holoEnsTestAcc(ii) = nanmean(holoCorrect(temp));
    end
else
    holoEnsTestAcc(1:length(ensIDs)) = -10;
    holoEnsCotuned(1:length(ensIDs)) = -10;
end



end

