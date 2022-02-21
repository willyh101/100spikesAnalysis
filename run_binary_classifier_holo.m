%% 
% Loads data from individual experiments to train a binary classifier
%
% The classifier is trained (and tested) on neuronal activity evoked by a
% holographic stimuli of 10 neurons
%
% loadList_all: function that loads the data 
% logistic_classifier: function that trains and tests the classifier
%
%
% Change cellsToInclude to perform the analysis on different groups of
% cells (e.g., default, notStimmed)
% Options for the cells to use: default, all, notStimmedOrOT, notStimmed, OT,
% stimmedAndOT, stimmed
%% Setup 
clear; close all; clc;

% addpath(genpath('100spikesAnalysis'))

restoredefaultpath;
folder = fileparts(which('run_binary_classifier_holo.m'));
addpath(genpath(folder));
rmpath(folder)

%% load data

loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/'; % where ever your files are
loadList_all = oriLoadList_GH('long');

%% Preallocation

% These are the angle pairings for the classifier
temp = [0:45:315]';
temp2 = [];
counter = 1;
for ii = 1:8
    for jj = 2:(8-ii+1)
        
        if abs(temp(1) - temp(jj))~=180
            temp2(counter,1) = temp(1);
            temp2(counter,2) = temp(jj);
        else
            temp2(counter,1) = nan;
            temp2(counter,2) = nan;
        end
        counter  = counter+1;
    end
    temp(1) = [];  
end
temp2(isnan(temp2(:,1)),:)=[];
ortho_pairings = temp2;

% ortho_pairings = mod([0:45:315; 90:45:405]',360);

holoTestAcc = nan(length(ortho_pairings),length(loadList_all));

%% Loop through the loadLists
for outer_loop = 1:length(loadList_all)
    
    %% Load the data (code from Will)
    loadList = loadList_all(outer_loop);
    [All,opts,outVars] = loadData_GH(loadList,loadPath);
    
    %% Get orientation tuning and OSI from dataset
    
    [All, outVars] = getTuningCurve(All, opts, outVars);
    [All, outVars] = calcOSI(All, outVars);
    [All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
    [All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically
    [All, outVars] = defineDistanceTypes(All, outVars);
    
    All.out.exp.ensMaxD =  outVars.ensMaxD;
    
    %% Run the classifier code
    
    %% Loop through the angle pairings
    fprintf('******************************\n')
    fprintf('Training and testing the classifier for different angle pairs for exp %d out of %d \n',outer_loop,length(loadList_all))
    for gg = 1:size(ortho_pairings,1)
        
        oriToUse = ortho_pairings(gg,:);
        % Only attempt if the visual stimulus was presented
        if sum(All.out.exp.ensPO(outVars.ensemblesToUse)==oriToUse(1))>0 ...
                && sum(All.out.exp.ensPO(outVars.ensemblesToUse)==oriToUse(2))>0
            
            % Options: default, all, notStimmedOrOT, notStimmed, OT,
            % stimmedAndOT, stimmed
            cellsToInclude = 'default';
            % Train and test the classifier
            [holoTestAcc(gg,outer_loop)]...
                = binary_classifier_holo(outVars,All,oriToUse,cellsToInclude);
           
            fprintf('Loop %d done\n',gg)
        else
            fprintf('Loop %d done: Not enough trials with these angles\n',gg)
        end 
    end
end


%% Plot the results
figure(1); clf; 
subplot(1,2,1); hold on;
boxplot(holoTestAcc(~isnan(holoTestAcc)))
plot([0.5 1.5], 50+0*[0.5 1.5],'k--')
set(gca,'fontsize',16)
ylabel('Testing Accuracy')
ylim([0 105])

subplot(1,2,2); hold on
ave_corr = nanmean(holoTestAcc(~isnan(holoTestAcc)));
err_corr = nanstd(holoTestAcc(~isnan(holoTestAcc)))/sqrt(sum(~isnan(holoTestAcc(:))));

errorbar(1,ave_corr,err_corr,'color','k','linewidth',1.5)
plot(1,ave_corr,'.','markersize',16,'color','k');
plot([0.5 1.5], 50+0*[0.5 1.5],'k--')
ylim([40 105])
set(gca,'fontsize',16)
ylabel('Testing Accuracy')
xticks([1])

sgtitle(sprintf(strcat(cellsToInclude, ' cells')),'fontsize',20,'FontWeight','bold')


