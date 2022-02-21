%%
% Trains (and tests) a binary classifier on neuronal activity evoked by
% holographic stimulations 
%%
function [holoTestAcc]=binary_classifier_holo(outVars,All,oriToUse,includedCells)

ensIDs = find((All.out.exp.ensPO==oriToUse(1) | All.out.exp.ensPO==oriToUse(2))  ...
    & outVars.ensemblesToUse==1);
ensOri = All.out.exp.ensPO(ensIDs);
ensIDs = [ensIDs(find(ensOri == oriToUse(1),1)) ensIDs(find(ensOri == oriToUse(2),1))];

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
    
    % trials to use
    trialsToUse_matrix(ii,:) = ismember(All.out.exp.stimID, u) &...
        All.out.exp.lowMotionTrials &...
        All.out.exp.lowRunTrials &...
        All.out.exp.stimSuccessTrial &...
        ismember(All.out.exp.visID, [0 1]);
    
    ensMaxD_trials(ii,:) = trialsToUse_matrix(ii,:)*All.out.exp.ensMaxD(ensID);
    orisTrialsMatrix(ii,:) = trialsToUse_matrix(ii,:)*All.out.exp.ensPO(ensID);
    
    % cells to use
    holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
    % tg = All.out.exp.rois{holo}; % this is cells in ensemble
    % fixed finding the cells in the ensembles?
    tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));
    
    offTargetRisks = All.out.anal.offTargetRisk(holo,:);
    ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
    visCells = All.out.anal.pVisR < 0.05; % visually responsive
    
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

% Store the zdfData from the right cells and trials
% cells x time x trials
% note: these are not trialwise baselined (but could be)
% traces = All.out.exp.dataToUse(cellsToUse, :, holoTrialsToUse);
tracesHolo = All.out.exp.zdfData(cellsToUse, :, holoTrialsToUse);

%% Optional pre-processing: zero-mean before stimulus across trials
baselineActHolo = squeeze(mean(mean(tracesHolo(:,1:6,:),2)));
for ii = 1: size(tracesHolo,3)
    tracesHolo(:,:,ii) = tracesHolo(:,:,ii)-baselineActHolo(ii);
end


%% Create the training and testing matrices

num_training_trials = min(floor(sum(ismember(holoOrisTrials,oriToUse(1)))/2),...
    floor(sum(ismember(holoOrisTrials,oriToUse(2)))/2));
num_testing_trials = num_training_trials;

if num_training_trials>7
    
    opts.recWinRange = [0.5 1.5];
    % timeframes for the visual stimulus
    recWinSec=opts.recWinRange+All.out.vis.visStart;
    
    % Use the activity when the visual stim was present
    winToUse = round(recWinSec*All.out.info.FR);
    % winToUse = [1 size(traces,2)]; % Alternative; use all data
    
    trace_ave = squeeze(mean(tracesHolo(:,winToUse(1):winToUse(2),:),2));
    training_matrix = [];
    trainingLabels = [];
    testing_matrix=[];
    testingLabels=[];
    for ii = 1:length(oriToUse)
        trainingTrials = find(holoOrisTrials==oriToUse(ii),num_training_trials);
        training_matrix = [training_matrix; trace_ave(:,trainingTrials)'];
        trainingLabels = [trainingLabels; oriToUse(ii)*ones(num_training_trials,1)];
        
        testingTrials = find(holoOrisTrials==oriToUse(ii),num_testing_trials,'last');
        testing_matrix = [testing_matrix; trace_ave(:,testingTrials)'];
        testingLabels = [testingLabels; oriToUse(ii)*ones(num_testing_trials,1)];
    end
    
    %% Train and test the classifier
    binary_classifier =fitclinear(training_matrix,trainingLabels,'Learner',...
        'svm');
    predictions = binary_classifier.predict(testing_matrix);
    
    % Save the testing accuracy for the classifier
    holoTestAcc = sum(predictions == testingLabels)/length(predictions)*100;
    
else
    holoTestAcc = nan;
end

end

