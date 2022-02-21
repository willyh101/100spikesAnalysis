%% 
% Loads data from individual experiments to train a binary classifier
%
% The classifier is trained (and tested) on neuronal activity evoked by a
% visual stimulus, and then tested on holographic stimulations
%
% loadList_all: function that loads the data 
% logistic_classifier: function that trains and tests the classifier
%%
%% setup junk
clear; close all; clc;

tic;

% addpath(genpath('100spikesAnalysis'))

restoredefaultpath;
folder = fileparts(which('run_binary_classifier.m')); 
addpath(genpath(folder));
rmpath(folder)

%% load data

% basePath = 'T:\Outfiles'; % where ever your files are
basePath = '/Users/gregoryhandy/Research_Local/outputdata1';

loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/';

% Smaller load list to debug the code
% loadList_all = {'210428_w32_2_outfile.mat'
%     '210903_I151_3_outfile.mat'
%     '210826_i151_3_outfile.mat'
%     '210914_i151_3_outfile.mat'
%     '210927_I154_2_outfile.mat'
%     '211019_I156_1_outfile.mat'};

% loadList_all = {'210428_w32_2_outfile.mat'
%     '210903_I151_3_outfile.mat'};


loadList_all = {'190418_I127_1_outfile.mat'
    '190506_I127_1_outfile.mat'
    '190513_I127_1_outfile.mat'
    '190515_I127_1_outfile.mat'
    '190525_I127_1_outfile.mat'
    '190529_I127_1_outfile.mat'
    '190604_I127_1_outfile.mat'
%     '191017_I136_1_outfile.mat' % something off with the # of time frames in tracesVis
%     '191126_I138_1_outfile.mat' %rate exp, but several valid; vis stim and holo stim off by 1 total frame
    '191206_I138_1_outfile.mat'
    '191217_mora_tre_outfile.mat'
    '200302_W14_1_outfile.mat'
    '200304_W14_1_outfile.mat'
    '200309_w21_1_outfile.mat'
    '200311_i139_2_outfile.mat'
%     '200728_i140_2_outfile.mat' % vis stim and holo stim off by 1 total frame (31 vs. 30)
    '200723_i140_2_outfile.mat'
    '200729_I140_2_outfile.mat'
    '200810_i140_2_outfile.mat'
    '200902_w18_3_outfile.mat'
    '200901_w19_1_outfile.mat'
    '201103_w29_1_outfile.mat'
%    '201112_w29_3_outfile.mat' % vis stim and holo stim off by 1 total frame
    '201116_w29_3_outfile.mat'
%    '201202_w29_3_outfile.mat' % vis stim and holo stim off by 1 total frame
    '210428_w32_2_outfile.mat' % first of optimizer
    '210902_I151_3_outfile.mat' %Many random holos at 30hz for SM Proj (strangely many failure..?)
    '210903_I151_3_outfile.mat'
    '210826_i151_3_outfile.mat'
    '210914_i151_3_outfile.mat'
    '210927_I154_2_outfile.mat'
    '211021_W40_2_outfile.mat' %rate Exp, but several valid
    '211019_I156_1_outfile.mat' %sepW1
    '211102_I158_1_outfile.mat' %sepW1
    '211108_I156_1_outfile.mat' %sepW1  one at a time
    '211019_I154_1_outfile.mat' %one at a time but some eligible
    };

overall_correct_all_outer = [];
testAcc_all_outer = [];
ensSizes = [];
class_trial_correct_all_outer = [];
class_trial_cotuned_all_outer = [];

%% Loop through the loadLists
for outer_loop = 1:length(loadList_all)

    %% Load the data (code from Will)
    loadList = loadList_all(outer_loop);
    [All,opts,outVars] = loadData_GH(loadList,loadPath);

    %% Orientation Tuning and OSI
    
    [All, outVars] = getTuningCurve(All, opts, outVars);
    [All, outVars] = calcOSI(All, outVars);
    [All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
    [All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically
    [All, outVars] = defineDistanceTypes(All, outVars);
    
    All.out.exp.ensMaxD =  outVars.ensMaxD;

%% Run the classifier code


% These are the angle pairings for the classifier
ortho_pairings = zeros(6,2);
ortho_pairings = [0 90];
for ii = 2:6
    ortho_pairings = [ortho_pairings; ortho_pairings(ii-1,:)+45];
end

ensIDs_all = [];
class_trial_correct_all = [];
class_trial_cotuned_all = [];
overall_correct_all = [];
pref_correct_all=[];
ortho_correct_all=[];
testAcc_all = [];
%% Loop through the angle pairings
fprintf('Training and testing the classifier for different angle pairs for exp %d out of %d \n',outer_loop,length(loadList_all)) 
for gg = 1:6
    
    oriToUse = ortho_pairings(gg,:);
    % Only attempt if the visual stimulus was presented
    if sum(All.out.exp.ensPO==oriToUse(1))>0 && sum(All.out.exp.ensPO==oriToUse(2))>0
        
        % Train and test the classifier
        [overall_correct,pref_correct,ortho_correct,class_trial_correct,...
            class_trial_cotuned,ensIDs,testAcc]...
            = logistic_classifier_v2(outVars,All,oriToUse);
        
        % Store the statistics
        overall_correct_all = [overall_correct_all,mean(overall_correct)];
        pref_correct_all = [pref_correct_all,mean(pref_correct)];
        ortho_correct_all = [ortho_correct_all,mean(ortho_correct)];
        
        testAcc_all = [testAcc_all,mean(testAcc(testAcc>70))];
        
        ensIDs_all = [ensIDs_all,ensIDs];
        class_trial_correct_all = [class_trial_correct_all,class_trial_correct];
        class_trial_cotuned_all = [class_trial_cotuned_all,class_trial_cotuned];
        fprintf('Loop %d done\n',gg)
    else
        fprintf('Loop %d done: Not enough trials with these angles\n',gg)
    end
    
    
end

% Store the overall statistics
overall_correct_all_outer = [overall_correct_all_outer, overall_correct_all];
testAcc_all_outer = [testAcc_all_outer,testAcc_all];

%% Plot the results with ensemble statistics
ensSizes = [ensSizes,All.out.exp.ensMaxD(ensIDs_all)];
class_trial_correct_all_outer = [class_trial_correct_all_outer, class_trial_correct_all];
class_trial_cotuned_all_outer = [class_trial_cotuned_all_outer,class_trial_cotuned_all];
% figure(32); hold on;
% scatter(All.out.exp.ensMaxD(ensIDs_all),class_trial_correct_all,[],class_trial_cotuned_all','filled')
% colormap(winter)
% colorbar
% caxis([0 1])
% set(gca,'fontsize',16)
% xlabel('Max Ensemble Distance')
% ylabel('Classifier accuracy')
end

%% Plot the results

figure(1243); clf; hold on;

for ii = 1:length(testAcc_all_outer)
    % Only plot the cases where the classifier is 'decent' at classifying
    % the visual stimuli
    if testAcc_all_outer(ii) >= 70
        plot([1 2],[testAcc_all_outer(ii) overall_correct_all_outer(ii)],'k-')
        plot([1 2],[testAcc_all_outer(ii) overall_correct_all_outer(ii)],'k.','markersize',16)
    end
end
set(gca,'fontsize',16)
xlim([0.8 2.2])
boxplot([testAcc_all_outer(testAcc_all_outer>70)';...
    overall_correct_all_outer(testAcc_all_outer>70)'],...
    [1*ones(length(testAcc_all_outer(testAcc_all_outer>70)),1);...
    2*ones(length(testAcc_all_outer(testAcc_all_outer>70)),1)])
ylim([0 105])
xticklabels({'Visual Stim','Holo Stim'})
ylabel('Testing Accuracy')


figure(124); clf;
subplot(1,2,1)
scatter(ensSizes,class_trial_correct_all_outer,[],class_trial_cotuned_all_outer','filled')
colormap(winter)
colorbar
caxis([0 1])
set(gca,'fontsize',16)
xlabel('Max Ensemble Distance')
ylabel('Classifier accuracy')


subplot(1,2,2)
boxplot([class_trial_correct_all_outer(ensSizes<350)'; ...
    class_trial_correct_all_outer(ensSizes>=350)'],...
    [1*ones(length(class_trial_correct_all_outer(ensSizes<350)),1);...
    2*ones(length(class_trial_correct_all_outer(ensSizes>=350)),1)])
set(gca,'fontsize',16)
xticks([1 2])
xticklabels({'Tight','Loose'})


figure(53); clf;
subplot(1,2,1)
scatter(ensSizes(class_trial_cotuned_all_outer==1),...
    class_trial_correct_all_outer(class_trial_cotuned_all_outer==1),'filled')
set(gca,'fontsize',16)
xlabel('Max Ensemble Distance')
ylabel('Classifier accuracy')
xlim([0 max(xlim)])

subplot(1,2,2)
boxplot([class_trial_correct_all_outer(ensSizes<350 & class_trial_cotuned_all_outer==1)';...
    class_trial_correct_all_outer(ensSizes>=350 & class_trial_cotuned_all_outer==1)'],...
    [1*ones(length(class_trial_correct_all_outer(ensSizes<350 & class_trial_cotuned_all_outer==1)),1);...
    2*ones(length(class_trial_correct_all_outer(ensSizes>=350 & class_trial_cotuned_all_outer==1)),1)])
set(gca,'fontsize',16)
xticks([1 2])
ylabel('Classifier accuracy')
xticklabels({'Tight','Loose'})


toc;