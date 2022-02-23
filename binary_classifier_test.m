clear; close all; clc;

%% Load in data from a single day
load('/Users/gregoryhandy/Research_Local/outputdata1/210903_I151_3_outfile.mat')

%% Parameters used to figure out which frames the vis. stim. is presented

opts.recWinRange = [0.5 1.5];
% timeframes for the visual stimulus 
recWinSec=opts.recWinRange+out.vis.visStart;
winToUse = round(recWinSec*out.info.FR);
bwinToUse = max(round([0 out.vis.visStart]*out.info.FR),[1 1]);

%% Find trials where the mouse is running and remove them
opts.runThreshold = 6 ; % trials with runspeed below this will be excluded
opts.runValPercent = 0.75; % percent of frames that need to be below run threshold

out.vis.lowRunTrials = ones(1,size(out.vis.runVal,1));
for trial_num = 1:size(out.vis.runVal,1)
    
    percent_run = length(find(out.vis.runVal(trial_num,1:winToUse(2))<opts.runThreshold))...
        /length(out.vis.runVal(trial_num,1:winToUse(2)));
    
    if percent_run<=opts.runValPercent
        out.vis.lowRunTrials(trial_num) = 0;
    end
end

%% Preallocate training matrix
num_neurons = size(out.vis.zdfData,1);
total_frames = size(out.vis.zdfData,2);

% Could adjust the start and and frame to be centered more on time frames
% when the vis. stim. is presented
start_frame = 1;
end_frame = total_frames;
num_frames = end_frame-start_frame+1;
num_trials = size(out.vis.zdfData,3);


num_training_trials = 5;
num_testing_trials = 7;
training_matrix = zeros(num_training_trials*2,num_neurons);
test_matrix = zeros(num_testing_trials*2,num_neurons);

%% Create the training matrix
% (total_frames x num of angles) x num_neurons

% visID: 2 is 0 degrees, 4 is 90 degrees, and 6 is 180
% Can adjust this to look at other angles as well
% visIDToTest = [2 4];
visIDToTest = [1 2];

vis_vals = zeros(length(visIDToTest),num_neurons);

for nn = 1:num_neurons
    for ii = 1:length(visIDToTest)
        wanted_trials = find(out.vis.visID==visIDToTest(ii) &...
            out.vis.lowMotionTrials &...
            out.vis.lowRunTrials);
        neuron_trials = squeeze(out.vis.zdfData(nn,1:num_frames,wanted_trials));
                 
        trial_ave = mean(neuron_trials,2);
        % Used for sorting by response during stim below
        vis_vals(ii,nn) = mean(trial_ave(winToUse(1):winToUse(2)));
        
        vis_ave = mean(neuron_trials(winToUse(1):winToUse(2),:));
        
        training_matrix(1+(ii-1)*num_training_trials:num_training_trials*ii,nn)...
            =vis_ave(1:num_training_trials);
        
        test_matrix(1+(ii-1)*num_testing_trials:num_testing_trials*ii,nn)...
            =vis_ave(1+num_training_trials:(num_training_trials+num_testing_trials));
    end
end

%% Create and test the classifier

num_to_use = [5 10 20 30 300];
per_corr = zeros(length(num_to_use),10);
for cc = 1:length(num_to_use)
    
    for inner_loop = 1:20
        winning_neurons = randperm(num_neurons,num_to_use(cc));
        
        adj_training_matrix=training_matrix(:,winning_neurons);
        adj_testing_matrix =test_matrix(:,winning_neurons);
        
        Y = zeros(num_training_trials*2,1);
        Y(1:num_training_trials) = 1;
        
        log_classifier =fitclinear(adj_training_matrix,Y,'Learner','logistic','Regularization','lasso');

        
        predictions = log_classifier.predict(adj_testing_matrix);
        correct = 0;
        for ii = 1:length(predictions)
            if predictions(ii) == 1 && ii <= num_testing_trials
                correct = correct + 1;
            elseif predictions(ii) == 0 && ii > num_testing_trials
                correct = correct + 1;
            end
        end
        per_corr(cc,inner_loop) = correct/length(predictions)*100;
    end
end

%%

figure(1); clf; hold on;
errorbar(mean(per_corr,2), std(per_corr,[],2)/sqrt(size(per_corr,2)))
plot(mean(per_corr,2),'.-','linewidth',1.5,'markersize',16,'color',[0 0.447 0.741])

xticks([1:length(num_to_use)])
xticklabels(num_to_use)
xlabel('Number of neurons used as predictors')
ylabel('Percent correct')
set(gca,'fontsize',16)
xlim([0.9 (length(num_to_use)+0.1)])

%%
% figure(3); clf; 
% subplot(1,2,1); hold on;
% imagesc(training_matrix(1:num_training_trials,winning_neurons)')
% caxis([0 1])
% xlim([0.5 (num_training_trials+0.5)])
% ylim([1 length(winning_neurons)])
% title('0 degree grating')
% set(gca,'fontsize',16)
% ylabel('Neuron number')
% xlabel('Trial number')
% 
% subplot(1,2,2); hold on;
% imagesc(training_matrix(num_training_trials+1:end,winning_neurons)')
% caxis([0 1])
% xlim([0.5 (num_training_trials+0.5)])
% ylim([1 length(winning_neurons)])
% title('90 degree grating')
% set(gca,'fontsize',16)
% yticks([])
% xlabel('Trial number')
% colorbar



% %% Perform PCA with built-in function
% 
% % Note: T*W' = training_matrix
% % So, T = training_matrix*W, since W'*W = I
% % W can be used to predict the trajectory in PC space of our stim trials
% [W, T, eigenvalues] = pca(training_matrix(:,pref_indices));
% 
% % Plot of percent variance explained
% figure(1); clf; hold on;
% plot(cumsum(eigenvalues/sum(eigenvalues)),'linewidth',1.5)
% set(gca,'fontsize',16)
% ylabel('Percenter variance explained')
% xlabel('Number of Components')
% 
% %% Plot the PCs for each stim angle
% 
% figure(2); clf; 
% subplot(1,3,1); hold on;
% maxVal = 0; minVal = 0;
% for ii = 1:length(visIDToTest)
%     maxVal = max(maxVal,max(T(1+num_frames*(ii-1):num_frames*ii,1)));
%     minVal = min(minVal,min(T(1+num_frames*(ii-1):num_frames*ii,1)));
%     plot(T(1+num_frames*(ii-1):num_frames*ii,1),'linewidth',1.5)
% end
% plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
% plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
% set(gca,'fontsize',16)
% xlabel('Time frames')
% ylabel('PC1')
% ylim([minVal maxVal])
% 
% subplot(1,3,2); hold on;
% maxVal = 0; minVal = 0;
% for ii = 1:length(visIDToTest)
%     maxVal = max(maxVal,max(T(1+num_frames*(ii-1):num_frames*ii,2)));
%     minVal = min(minVal,min(T(1+num_frames*(ii-1):num_frames*ii,2)));
%     plot(T(1+num_frames*(ii-1):num_frames*ii,2),'linewidth',1.5)
% end
% set(gca,'fontsize',16)
% xlabel('Time')
% ylabel('PC2')
% plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
% plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
% 
% subplot(1,3,3); hold on;
% for ii = 1:length(visIDToTest)
%     h_leg(ii) = plot(T(1+num_frames*(ii-1):num_frames*ii,1),T(1+num_frames*(ii-1):num_frames*ii,2),'linewidth',1.5);
%     plot(T(1+num_frames*(ii-1),1),T(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
% end
% 
% set(gca,'fontsize',16)
% xlabel('PC1')
% ylabel('PC2')
% 
% legend(h_leg,{'0 degree','90 degree','180 degree'})
% 
