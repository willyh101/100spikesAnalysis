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

training_matrix = zeros(num_frames*2,num_neurons);

%% Create the training matrix
% (total_frames x num of angles) x num_neurons

% visID: 2 is 0 degrees, 4 is 90 degrees, and 6 is 180
% Can adjust this to look at other angles as well
visIDToTest = [2 4 6];

vis_vals = zeros(length(visIDToTest),num_neurons);

for nn = 1:num_neurons
    for ii = 1:length(visIDToTest)
        wanted_trials = find(out.vis.visID==visIDToTest(ii) &...
            out.vis.lowMotionTrials &...
            out.vis.lowRunTrials);
        neuron_trials = squeeze(out.vis.zdfData(nn,1:num_frames,wanted_trials));
        
        % Alternative: z-score the allData matrix by trial
        % Same story with either method
%         neuron_trials = squeeze(out.vis.allData(nn,start_frame:end_frame,wanted_trials)); 
%         for ll = 1:size(neuron_trials,2)      
%             ave_est=mean(neuron_trials(:,ll));
%             sigma_est=std(neuron_trials(:,ll));
%             neuron_trials(:,ll) = (neuron_trials(:,ll)-ave_est)/sigma_est;
%         end
         
        trial_ave = mean(neuron_trials,2);
        
        % Used for sorting by response during stim below
        vis_vals(ii,nn) = mean(trial_ave(winToUse(1):winToUse(2)));
        training_matrix(1+(ii-1)*num_frames:num_frames*ii,nn)=trial_ave;
    end
end

%% Keep just the 200 strongest responding neurons to 0 and 90 degrees


[~,I] = sort(vis_vals(1,:),'descend');
pref0_indices = I(1:200);
[~,I] = sort(vis_vals(2,:),'descend');
pref90_indices = I(1:200);
[~,I] = sort(vis_vals(3,:),'descend');
pref180_indices = I(1:200);

pref_indices = unique([pref0_indices pref90_indices pref180_indices],'stable');

figure(3); clf; 
subplot(1,3,1); hold on;
imagesc(training_matrix(1:num_frames,pref_indices)')
caxis([0 1])
xlim([1 total_frames])
ylim([1 length(pref_indices)])
title('0 degree grating')
set(gca,'fontsize',16)
ylabel('Neuron number')
xlabel('Time frames')

subplot(1,3,2); hold on;
imagesc(training_matrix(1+num_frames:num_frames*2,pref_indices)')
caxis([0 1])
xlim([1 total_frames])
ylim([1 length(pref_indices)])
title('90 degree grating')
set(gca,'fontsize',16)
yticks([])
xlabel('Time frames')

subplot(1,3,3); hold on;
imagesc(training_matrix(1+num_frames*2:num_frames*3,pref_indices)')
colorbar
caxis([0 1])
yticks([])
xlim([1 total_frames])
ylim([1 length(pref_indices)])
title('180 degree grating')
xlabel('Time frames')
set(gca,'fontsize',16)


%% Perform PCA with built-in function

% Note: T*W' = training_matrix - mean(training_matrix)
% So, T = training_matrix*W, since W'*W = I
% W can be used to predict the trajectory in PC space of our stim trials
[W, T, eigenvalues] = pca(training_matrix(1:2*num_frames,pref_indices));

% Plot of percent variance explained
figure(1); clf; hold on;
plot(cumsum(eigenvalues/sum(eigenvalues)),'linewidth',1.5)
set(gca,'fontsize',16)
ylabel('Percenter variance explained')
xlabel('Number of Components')

%%

test_matrix = training_matrix(num_frames*2+1:end,pref_indices)-...
    mean(training_matrix(num_frames*2+1:end,pref_indices));

projection_180 = test_matrix*W;

%% Plot the PCs for each stim angle

figure(2); clf; 
subplot(1,3,1); hold on;
maxVal = 0; minVal = 0;
for ii = 1:length(visIDToTest)-1
    maxVal = max(maxVal,max(T(1+num_frames*(ii-1):num_frames*ii,1)));
    minVal = min(minVal,min(T(1+num_frames*(ii-1):num_frames*ii,1)));
    plot(T(1+num_frames*(ii-1):num_frames*ii,1),'linewidth',1.5)
end
plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')

plot(projection_180(:,1),'linewidth',1.5)
set(gca,'fontsize',16)
xlabel('Time frames')
ylabel('PC1')
ylim([minVal maxVal])

subplot(1,3,2); hold on;
maxVal = 0; minVal = 0;
for ii = 1:length(visIDToTest)-1
    maxVal = max(maxVal,max(T(1+num_frames*(ii-1):num_frames*ii,2)));
    minVal = min(minVal,min(T(1+num_frames*(ii-1):num_frames*ii,2)));
    plot(T(1+num_frames*(ii-1):num_frames*ii,2),'linewidth',1.5)
end

plot(projection_180(:,2),'linewidth',1.5)
set(gca,'fontsize',16)
xlabel('Time')
ylabel('PC2')
plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')

subplot(1,3,3); hold on;
for ii = 1:length(visIDToTest)-1
    h_leg(ii) = plot(T(1+num_frames*(ii-1):num_frames*ii,1),T(1+num_frames*(ii-1):num_frames*ii,2),'linewidth',1.5);
    plot(T(1+num_frames*(ii-1),1),T(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
end

h_leg(3) = plot(projection_180(:,1),projection_180(:,2),'linewidth',1.5);
plot(projection_180(1,1),projection_180(1,2),'k.','markersize',16,'linewidth',1.5)
set(gca,'fontsize',16)
xlabel('PC1')
ylabel('PC2')

legend(h_leg,{'0 degree','90 degree','180 degree'})

