clear; close all; clc;

% Adds all subfolders to the path
restoredefaultpath;
folder = fileparts(which('run_PCA_compare_holo.m')); 
addpath(genpath(folder));
rmpath(folder);


%% Load data
% Options: default, all, stimmedAndOT, stimmed, notStimmed
includedCells = 'all';

loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/';
% Load up single experiment
loadList_all = {'210903_I151_3_outfile.mat'};
% loadList_all = {'210927_I154_2_outfile.mat'};
% Load the data (code from Will)
loadList = loadList_all(1);
[All,opts,outVars] = loadData_GH(loadList,loadPath);

%% Orientation Tuning and OSI
[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
[All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically
[All, outVars] = defineDistanceTypes(All, outVars);

All.out.exp.ensMaxD =  outVars.ensMaxD;

%% Choose ensembles
oriToUse = [180 270];
% oriToUse = [0 90];

ensIDs  = zeros(1,2);
ensOri = zeros(1,2);
for ii = 1:2
    temp = find((All.out.exp.ensPO==oriToUse(ii)) & outVars.ensemblesToUse==1);
    if ~isempty(temp)
        ensIDs(ii) = temp(ceil(rand()*length(temp)));
        ensOri(ii) = All.out.exp.ensPO(ensIDs(ii));
    else
        error('No ensembles with the requested orientation')
    end
end

%% Timeframes for the visual stimulus 
opts.recWinRange = [0 1];
recWinSec=opts.recWinRange+All.out.vis.visStart;
winToUse = round(recWinSec*All.out.info.FR);
num_frames = 24;

%% Preallocation

numTotalTrials = length(All.out.exp.lowMotionTrials);
numTotalCells = length(All.out.exp.allDepth);

trialsToUse_matrix=zeros(length(ensIDs), numTotalTrials);
orisTrialsMatrix = zeros(length(ensIDs), numTotalTrials);
cellsToUse_matrix = zeros(length(ensIDs), numTotalCells);
tg_matrix = zeros(length(ensIDs), numTotalCells);

%%
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
    
    orisTrialsMatrix(ii,:) = ensOri(ii)*trialsToUse_matrix(ii,:);
    % cells to use
    holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
    tg = All.out.exp.rois{holo}; % this is cells in ensemble
    tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));

    offTargetRisks = All.out.anal.offTargetRisk(holo,:);
    ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
    visCells = All.out.anal.pVisR < 0.05; % visually responsive
    
    if strcmp(includedCells,'default')
        cellsToUse_matrix(ii,:) = ~offTargetRisks &...
            ~ROIinArtifact' &...
            visCells;
    elseif strcmp(includedCells,'all')
        cellsToUse_matrix(ii,:) = ones(size(offTargetRisks));
    elseif strcmp(includedCells,'notStimmed')
        cellsToUse_matrix(ii,:) = ~ismember(1:numel(offTargetRisks), tg_adj);
    elseif strcmp(includedCells,'stimmedAndOT')
        cellsToUse_matrix(ii,:) = offTargetRisks;
    elseif strcmp(includedCells,'stimmed')
        cellsToUse_matrix(ii,:) = ismember(1:numel(offTargetRisks), tg_adj);
    else
        error('Unknown cell selection');
    end
end

if strcmp(includedCells,'stimmedAndOT') || strcmp(includedCells,'stimmed')
    cellsToUse = sum(cellsToUse_matrix,1)>0;
else % only include cells that were never stimmed/off-targets
    cellsToUse = logical(prod(cellsToUse_matrix));
end
holoTrialsToUse = logical(sum(trialsToUse_matrix,1));
holoOrisTrials = sum(orisTrialsMatrix,1);
holoOrisTrials = holoOrisTrials(holoTrialsToUse);

tracesHolo = All.out.exp.zdfData(cellsToUse, :, holoTrialsToUse);

%% Split the trials into two (training/testing), and average across them
holoTestMatrix = zeros(num_frames*2,size(tracesHolo,1));
holoMatrix = zeros(num_frames*2,size(tracesHolo,1));
% average across trials
for ii = 1:length(ensIDs)
    
    tempArray = find(holoOrisTrials==oriToUse(ii));
    totalTrials = length(tempArray);
    if totalTrials>10
        firstHalf = floor(totalTrials/2);
        tempArray = tempArray(randperm(totalTrials));
        
        trialAveHolo = mean(tracesHolo(:,:,tempArray(1:firstHalf)),3);
        holoMatrix(1+(ii-1)*num_frames:ii*num_frames,:) = trialAveHolo';
        
        trialAveHolo = mean(tracesHolo(:,:,tempArray(1+firstHalf:end)),3);
        holoTestMatrix(1+(ii-1)*num_frames:ii*num_frames,:) = trialAveHolo';
    else
        error('Probably not enough trials to proceed')
    end
end

%% Optional pre-processing (zero-mean before stimulus across trials) 
for ii = 1:size(holoMatrix,2)
    for jj = 1:length(ensIDs)
        holoMatrix(1+num_frames*(jj-1):num_frames*jj,ii) ...
            = holoMatrix(1+num_frames*(jj-1):num_frames*jj,ii)-...
            mean(holoMatrix(1+num_frames*(jj-1):num_frames*(jj-1)+6,ii));
        
        holoTestMatrix(1+num_frames*(jj-1):num_frames*jj,ii) ...
            = holoTestMatrix(1+num_frames*(jj-1):num_frames*jj,ii)-...
            mean(holoTestMatrix(1+num_frames*(jj-1):num_frames*(jj-1)+6,ii));
    end
end

%%
holoAve = mean(holoMatrix(winToUse(1):winToUse(2),:));
[val,indices] = sort(holoAve,'descend');
%%
colorscheme = [0    0.4470    0.7410; 0.8500    0.3250    0.0980];

figure(11); clf; hold on;
aveResponse = mean(holoMatrix,2);
aveTestResponse = mean(holoTestMatrix,2);
for ii = 1:2
    g_final(ii)=plot(aveResponse(1+(ii-1)*num_frames:ii*num_frames),'color',colorscheme(ii,:),'linewidth',1.5);
    h_final(ii)=plot(aveTestResponse(1+(ii-1)*num_frames:ii*num_frames),'--','color',colorscheme(ii,:),'linewidth',1.5);
end
minVal = min(min(min(aveResponse)),min(min(aveTestResponse))); 
maxVal = max(max(max(aveResponse)),max(max(aveTestResponse)));
plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
ylim([minVal maxVal])
plot([1 num_frames],0*[1 num_frames],'k--')

set(gca,'fontsize',16)
xlabel('Time Frames')
ylabel('Normalized Response')
legend([g_final h_final],{strcat(sprintf('%0.0f',oriToUse(1)),'\circ vis'),...
    strcat(sprintf('%0.0f',oriToUse(2)),'\circ vis'),...
    strcat(sprintf('%0.0f',oriToUse(1)),'\circ vis'),...
    strcat(sprintf('%0.0f',oriToUse(2)),'\circ vis')})

%%

figure(2); clf; 
subplot(2,2,1); hold on;
imagesc(holoMatrix(1:num_frames,indices)')
caxis([0 1])
xlim([1 size(holoMatrix,1)/2])
ylim([1 size(holoMatrix,2)])
title(sprintf('%0.0f degree grating',oriToUse(1)))
set(gca,'fontsize',16)
ylabel('Neuron number')
xlabel('Time frames')
colorbar()

subplot(2,2,2); hold on;
imagesc(holoMatrix(1+num_frames:num_frames*2,indices)')
caxis([0 1])
xlim([1 size(holoMatrix,1)/2])
ylim([1 size(holoMatrix,2)])
title(sprintf('%0.0f degree grating',oriToUse(2)))
set(gca,'fontsize',16)
ylabel('Neuron number')
xlabel('Time frames')
colorbar()

subplot(2,2,3); hold on;
imagesc(holoTestMatrix(1:num_frames,indices)')
caxis([0 1])
xlim([1 size(holoMatrix,1)/2])
ylim([1 size(holoMatrix,2)])
title(sprintf('%0.0f degree grating',oriToUse(1)))
set(gca,'fontsize',16)
ylabel('Neuron number')
xlabel('Time frames')
colorbar()

subplot(2,2,4); hold on;
imagesc(holoTestMatrix(1+num_frames:num_frames*2,indices)')
caxis([0 1])
xlim([1 size(holoMatrix,1)/2])
ylim([1 size(holoMatrix,2)])
title(sprintf('%0.0f degree grating',oriToUse(2)))
set(gca,'fontsize',16)
ylabel('Neuron number')
xlabel('Time frames')
colorbar()



%% Perform PCA on half the trials and project the other half

[W_holo, T_holo, eigenvalues_holo] = pca(holoTestMatrix);

holoProjection = (holoMatrix-mean(holoMatrix))*W_holo;

figure(1); clf; hold on;
plot(cumsum(eigenvalues_holo/sum(eigenvalues_holo)),'linewidth',1.5)
set(gca,'fontsize',16)
ylabel('Percenter variance explained')
xlabel('Number of Components')
title('Holo PCA')

participation_ratio = sum(eigenvalues_holo)/sum(eigenvalues_holo.^2)


%%

figure(4); clf; 
subplot(1,3,1); hold on;
maxVal = 0; minVal = 0;
for ii = 1:2
    maxVal = max(maxVal,max(T_holo(1+num_frames*(ii-1):num_frames*ii,1)));
    minVal = min(minVal,min(T_holo(1+num_frames*(ii-1):num_frames*ii,1)));
    plot(T_holo(1+num_frames*(ii-1):num_frames*ii,1),'-','linewidth',1.5,'color',colorscheme(ii,:))
    plot(holoProjection(1+num_frames*(ii-1):num_frames*ii,1),'--','linewidth',1.5,'color',colorscheme(ii,:))
end
plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
set(gca,'fontsize',16)
xlabel('Time frames')
ylabel('PC1')
ylim([minVal maxVal])

subplot(1,3,2); hold on;
maxVal = 0; minVal = 0;
for ii = 1:2
    maxVal = max(maxVal,max(T_holo(1+num_frames*(ii-1):num_frames*ii,2)));
    minVal = min(minVal,min(T_holo(1+num_frames*(ii-1):num_frames*ii,2)));
    plot(T_holo(1+num_frames*(ii-1):num_frames*ii,2),'-','linewidth',1.5,...
        'color',colorscheme(ii,:))
    
    plot(holoProjection(1+num_frames*(ii-1):num_frames*ii,2),'--','linewidth',1.5,...
        'color',colorscheme(ii,:))
end
set(gca,'fontsize',16)
xlabel('Time')
ylabel('PC2')
plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
ylim([minVal maxVal])

subplot(1,3,3); hold on;
for ii = 1:2
    g_leg(ii) = plot(T_holo(1+num_frames*(ii-1):num_frames*ii,1),T_holo(1+num_frames*(ii-1):num_frames*ii,2),...
        '-','linewidth',1.5,'color',colorscheme(ii,:));
    plot(T_holo(1+num_frames*(ii-1),1),T_holo(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
    
    h_leg(ii) = plot(holoProjection(1+num_frames*(ii-1):num_frames*ii,1),holoProjection(1+num_frames*(ii-1):num_frames*ii,2),...
        '--','linewidth',1.5,'color',colorscheme(ii,:));
    plot(holoProjection(1+num_frames*(ii-1),1),holoProjection(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
    
end

set(gca,'fontsize',16)
xlabel('PC1')
ylabel('PC2')

legend([g_leg h_leg],{strcat(sprintf('%0.0f',oriToUse(1)),'\circ holo'),...
    strcat(sprintf('%0.0f',oriToUse(2)),'\circ holo'),...
    strcat(sprintf('%0.0f',oriToUse(1)),'\circ holo'),...
    strcat(sprintf('%0.0f',oriToUse(2)),'\circ holo')})




