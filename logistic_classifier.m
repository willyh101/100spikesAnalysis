%%
% Trains (and tests) a binary classifier on neuronal activity evoked by a
% visual stimulus, and then attempts to classify the ensemble orientation preference
% based on the activity evoked from holographic stimulation
%%
function [overall_correct,pref_correct,ortho_correct,class_trial_correct,...
    class_trial_cotuned,ensIDs,testAcc]=logistic_classifier(outVars,All,oriToUse)

ensIDs = [find(All.out.exp.ensPO==oriToUse(1)) find(All.out.exp.ensPO==oriToUse(2))];
ensOri = [oriToUse(1)*ones(size(find(All.out.exp.ensPO==oriToUse(1)),2),1);...
    oriToUse(2)*ones(size(find(All.out.exp.ensPO==oriToUse(2)),2),1)];

ensCoTuned = All.out.exp.meanEnsOSI(ensIDs)>0.2 & All.out.exp.ensOSI(ensIDs)>0.2;

%%

trialsToUse_matrix=[];
orisTrialsMatrix = [];
cellsToUse_matrix = [];
ensMaxD_trials = [];
ensCoTuned_trials = [];
for ii = 1:length(ensIDs)
    ens = ensIDs(ii);
    
    s = outVars.ensHNumber(ens);
    us = unique(All.out.exp.stimID);
    u = us(s);
    
    % trials to use
    trialsToUse_matrix(ii,:) = ismember(All.out.exp.stimID, u) &...
        All.out.exp.lowMotionTrials &...
        All.out.exp.lowRunTrials &...
        All.out.exp.stimSuccessTrial &...
        ismember(All.out.exp.visID, [0 1]);
    
    ensMaxD_trials(ii,:) = trialsToUse_matrix(ii,:)*All.out.exp.ensMaxD(ensIDs(ii));
    ensCoTuned_trials(ii,:) = trialsToUse_matrix(ii,:)*ensCoTuned(ii);
    
    orisTrialsMatrix(ii,:) = ensOri(ii)*trialsToUse_matrix(ii,:);
    % cells to use
    holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
    tg = All.out.exp.rois{holo}; % this is cells in ensemble
    offTargetRisks = All.out.anal.offTargetRisk(holo,:);
    ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
    visCells = All.out.anal.pVisR < 0.05; % visually responsive
    orisUsed = [nan 0:45:315];
    prefOris = arrayfun(@(x) orisUsed(x), All.out.anal.prefOri);
    
    % comment out parts if you want to change subsets of cells
    % note ~offTargetRisks excludes the targeted cells, if you want to keep the
    % off targets and exclude the targets (for this specific hologram) you can
    % use: ~ismember(1:numel(offTargetRisks), tg)
    cellsToUse_matrix(ii,:) = ~offTargetRisks &...
        ~ROIinArtifact' &...
        visCells;
%     cellsToUse_matrix(ii,:) = ~offTargetRisks;

    
%     disp('....')
%     disp(['Ensemble is tuned for ' num2str(All.out.exp.ensPO(ens))])
%     disp(['Ensemble OSI: ' num2str(All.out.exp.ensOSI(ens))])
%     disp(['Mean OSI: ' num2str(All.out.exp.meanEnsOSI(ens))])
%     disp(['Ensemble has a stim score of ' num2str(All.out.exp.ensStimScore(ens))])
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

% disp(['Total of ' num2str(sum(holoTrialsToUse)) ' trials used.'])
% disp(['Total of ' num2str(sum(cellsToUse)) ' cells used.'])

% here are you traces from the experiment
% cells x time x trials
% note: these are not trialwise baselined (but could be)
% traces = All.out.exp.dataToUse(cellsToUse, :, holoTrialsToUse);
traces = All.out.exp.zdfData(cellsToUse, :, holoTrialsToUse);

ensMaxD_trials = sum(ensMaxD_trials(:,holoTrialsToUse));
ensCoTuned_trials = sum(ensCoTuned_trials(:,holoTrialsToUse));

%% vis epoch
% if you want the responses from the vis epoch, for the same cells, run
% this code

% you can change to out.vis.dfData if it exists
visData = All.out.vis.zdfData;

% trials to use
vs = All.out.vis.visID;
vs(vs==0) = []; % Band-aid on an unknown issue; vs should be between 1 and 9
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

%% Create the training and testing matrices 

num_training_trials = min(floor(sum(ismember(oris, orisToUseVis(1)))/2),...
    floor(sum(ismember(oris, orisToUseVis(2)))/2));
num_testing_trials = num_training_trials;
num_neurons = size(tracesVis,1);

opts.recWinRange = [0.5 1.5];
% timeframes for the visual stimulus 
recWinSec=opts.recWinRange+All.out.vis.visStart;

% Use the activity when the visual stim was present
winToUse = round(recWinSec*All.out.info.FR);

% winToUse = [1 size(tracesVis,2)]; % Alternative; use all data
% bwinToUse = max(round([0 All.out.vis.visStart]*All.out.info.FR),[1 1]);
size(tracesVis)
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

%% Train 100 classifiers to assess success of classifying holo stim tuning
num_trainings = 100;
trial_correct = [];
for training_loop = 1:num_trainings

    %% Fit the logistic classifier on the training set of visual stims
    log_classifier =fitclinear(training_matrix,trainingLabels,'Learner',...
        'logistic','Regularization','lasso');
    
    %% Predict the testing set with classifier
    predictions = log_classifier.predict(testing_matrix);
    correct = 0;
    for ii = 1:length(predictions)
        if predictions(ii) == orisToUseVis(1) && ii <= num_testing_trials
            correct = correct + 1;
        elseif predictions(ii) == orisToUseVis(2) && ii > num_testing_trials
            correct = correct + 1;
        end
    end
    
    % Save the testing accuracy for the classifier
    testAcc(training_loop) = round(correct/length(predictions)*100);
    
    %% Predict holographic tuning with classifier
    holo_ave = squeeze(mean(traces(:,winToUse(1):winToUse(2),:),2))';
    holo_Predictions = log_classifier.predict(holo_ave);
    
    for ii = 1:length(holo_Predictions)
        if holo_Predictions(ii) == holoOrisTrials(ii)
            trial_correct(ii,training_loop) = 1;
        else
            trial_correct(ii,training_loop) = 0;
        end
    end
end
%%
overall_correct = sum(trial_correct)/size(trial_correct,1)*100;
pref_correct = sum(trial_correct(holoOrisTrials==oriToUse(1),:))/sum(holoOrisTrials==oriToUse(1))*100;
ortho_correct = sum(trial_correct(holoOrisTrials==oriToUse(2),:))/sum(holoOrisTrials==oriToUse(2))*100;

%%
class_trial_correct = [];
class_trial_cotuned = [];
for ii =1:length(ensIDs)
    temp = (ensMaxD_trials == All.out.exp.ensMaxD(ensIDs(ii)));
    
    if All.out.exp.meanEnsOSI(ensIDs(ii))>0.5 && All.out.exp.ensOSI(ensIDs(ii))>0.5
        class_trial_cotuned(ii) = 1;
    else
        class_trial_cotuned(ii) = 0;
    end
    class_trial_correct(ii) = mean(sum(trial_correct(temp,:),2)/num_trainings);
end

end

