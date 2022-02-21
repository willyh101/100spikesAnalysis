%% 
% Loads data from individual experiments to train a binary classifier
%
% The classifier is trained (and tested) on neuronal activity evoked by a
% visual stimulus, and then tested on holographic stimulations
%
% loadList_all: function that loads the data 
% logistic_classifier: function that trains and tests the classifier
%
% Change cellsToInclude to perform the analysis on different groups of
% cells (e.g., default, notStimmed)
%%
clear; close all; clc;

% Adds all subfolders to the path
restoredefaultpath;
folder = fileparts(which('run_binary_classifier_vis.m')); 
addpath(genpath(folder));
rmpath(folder);

%% load data

loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/'; % where ever your files are
loadList_all = oriLoadList_GH('long');

%%

ensSizes = [];
ensSizes_v2 = [];
holoEnsTestAccAll = [];
holoEnsCotunedAll = [];

%% Preallocation
% These are the angle pairings for the classifier
ortho_pairings = mod([0:45:315; 90:45:405]',360);

holoTestAcc = nan(length(ortho_pairings),length(loadList_all));
visTestAcc = nan(length(ortho_pairings),length(loadList_all));

%% Loop through the loadLists
for outer_loop = 1:length(loadList_all)
    
    %% Load the data from a specific experiment (code from Will)
    loadList = loadList_all(outer_loop);
    [All,opts,outVars] = loadData_GH(loadList,loadPath);
        
    
    %% Get orientation tuning and OSI from dataset
    
    [All, outVars] = getTuningCurve(All, opts, outVars);
    [All, outVars] = calcOSI(All, outVars);
    [All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
    [All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically
    [All, outVars] = defineDistanceTypes(All, outVars);
    All.out.exp.ensMaxD =  outVars.ensMaxD;
    
    %% Run the classifier code for all angle pairings
    fprintf('******************************\n')
    fprintf('Training and testing the classifier for different angle pairs for exp %d out of %d \n',outer_loop,length(loadList_all))
    for gg = 1:size(ortho_pairings,1)
        
        oriToUse = ortho_pairings(gg,:);
        % Only attempt if both visual stimuli were presented
        if sum(All.out.exp.ensPO(outVars.ensemblesToUse)==oriToUse(1))>0 ...
                && sum(All.out.exp.ensPO(outVars.ensemblesToUse)==oriToUse(2))>0
            
            % Options: default, all, notStimmedOrOT, notStimmed, OT,
            % stimmedAndOT, stimmed
            cellsToInclude = 'default';
            % Train and test the classifier
            [holoTestAcc(gg,outer_loop),visTestAcc(gg,outer_loop),...
                holoEnsTestAcc,holoEnsCotuned,ensIDs]=...
                binary_classifier_vis(outVars,All,oriToUse,'svm','',cellsToInclude);
            % Store the statistics            
            ensSizes = [ensSizes,All.out.exp.ensMaxD(ensIDs)];
            holoEnsTestAccAll = [holoEnsTestAccAll,holoEnsTestAcc];
            holoEnsCotunedAll = [holoEnsCotunedAll,holoEnsCotuned];
        end
    end  
end
%%
holoTestAcc = holoTestAcc(~isnan(holoTestAcc(:)));
visTestAcc=visTestAcc(~isnan(visTestAcc(:)));

holoEnsTestAccAll = holoEnsTestAccAll(~isnan(holoEnsTestAccAll));
holoEnsCotunedAll = holoEnsCotunedAll(~isnan(holoEnsCotunedAll));

%% Plot the results

figure(1243); clf; 

subplot(1,2,1);hold on;

for ii = 1:length(visTestAcc)
    % Only plot the cases where the classifier is 'decent' at classifying
    % the visual stimuli
    if visTestAcc(ii) >= 70
        plot([1 2],[visTestAcc(ii) holoTestAcc(ii)],'k-')
        plot([1 2],[visTestAcc(ii) holoTestAcc(ii)],'k.','markersize',16)
    end
end
set(gca,'fontsize',16)

boxplot([visTestAcc(visTestAcc>70);holoTestAcc(visTestAcc>70)],...
    [1*ones(length(visTestAcc(visTestAcc>70)),1);...
    2*ones(length(visTestAcc(visTestAcc>70)),1)])
ylim([0 105])
xticklabels({'Visual Stim','Holo Stim'})
ylabel('Testing Accuracy')

plot([0.5 2.5], 50+0*[0.8 2.2],'k--')
xlim([0.5 2.5])

box off


% figure(124); clf;
% subplot(1,2,1)
% scatter(ensSizes(holoEnsTestAccAll>=0),holoEnsTestAccAll(holoEnsTestAccAll>=0),'filled')
% set(gca,'fontsize',16)
% xlabel('Max Ensemble Distance')
% ylabel('Classifier accuracy')
% title('All Ensembles')
% 
% % class_trial_cotuned_all_outer =class_trial_cotuned_all_outer(~isnan(class_trial_cotuned_all_outer));
% 
% subplot(1,2,2)
% scatter(ensSizes(holoEnsCotunedAll==1),...
%     holoEnsTestAccAll(holoEnsCotunedAll==1),'filled')
% set(gca,'fontsize',16)
% xlabel('Max Ensemble Distance')
% ylabel('Classifier accuracy')
% xlim([0 max(xlim)])
% title('Co-tuned Ensembles')


% figure(54); clf;
subplot(1,2,2); hold on;

ave_vis = mean(visTestAcc(visTestAcc>70));
error_vis = std(visTestAcc(visTestAcc>70))/sqrt(length(visTestAcc));

errorbar(0,ave_vis,error_vis,'color','k','linewidth',1.5)
h = plot(0,ave_vis,'.','markersize',16,'color','k');


ensBins = [0 350 inf];
cotuned = [1 0];

color_scheme = [0    0.4470    0.7410; 0.8500    0.3250    0.0980];


gridSpacing = [0 0.3 0.6 1 1.3];
counter = 2;
for ii = 1:2
    for jj = 1:2
        trial_curr = holoEnsTestAccAll(ensSizes>ensBins(jj) &...
            ensSizes<ensBins(jj+1) & holoEnsCotunedAll==cotuned(ii));
        ave_corr = nanmean(trial_curr)*100;
        err_corr = nanstd(trial_curr)/sqrt(sum(~isnan(trial_curr)))*100;
        
        errorbar(gridSpacing(counter),ave_corr,err_corr,'color',color_scheme(ii,:),'linewidth',1.5)
        g(ii) = plot(gridSpacing(counter),ave_corr,'.','markersize',16,'color',color_scheme(ii,:));

        counter = counter+1;
    end
end

ylim([30 105])
xlim([-0.25 1.5])
set(gca,'fontsize',16)
xticks(gridSpacing)
xticklabels({'','Tight','Loose','Tight','Loose'})

plot([-0.25 1.75],[50 50],'k--')
legend([h g], {'Vis','Co-tuned','Mixed tuned'})
ylabel('Percent correct (classifier)')

sgtitle(sprintf(strcat(cellsToInclude, ' cells')),'fontsize',20,'FontWeight','bold')

