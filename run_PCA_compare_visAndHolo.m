clear; close all; clc;

%% Load data
% Options: default, all, stimmedAndOT, stimmed
includedCells = 'stimmed';

loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/';
% Load up single experiment
loadList_all = {'210903_I151_3_outfile.mat'};
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

%% Now get the data from vis epochs

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
tracesVis = visData(cellsToUse, :, trialsToUse);
trialAngle = oris(trialsToUse);

%% Average the time courses across trials
holoMatrix = zeros(num_frames*length(oriToUse),size(tracesVis,1));
visMatrix = zeros(num_frames*length(oriToUse),size(tracesVis,1));
for ii = 1:length(oriToUse)
    trialAveHolo = mean(tracesHolo(:,:,holoOrisTrials==oriToUse(ii)),3);
    holoMatrix(1+(ii-1)*num_frames:ii*num_frames,:) = trialAveHolo';
    
    trialAveVis = mean(tracesVis(:,:,trialAngle==oriToUse(ii)),3);
    visMatrix(1+(ii-1)*num_frames:ii*num_frames,:) = trialAveVis';
end

%% Optional pre-processing (zero-mean before stimulus across trials)

for ii = 1:size(holoMatrix,2)
    for jj = 1:length(ensIDs)
        holoMatrix(1+num_frames*(jj-1):num_frames*jj,ii) ...
            = holoMatrix(1+num_frames*(jj-1):num_frames*jj,ii)-...
            mean(holoMatrix(1+num_frames*(jj-1):num_frames*(jj-1)+6,ii));
    end
end

for ii = 1:size(visMatrix,2)
    visMatrix(1:num_frames,ii) = visMatrix(1:num_frames,ii)-mean(visMatrix(1:6,ii));
    visMatrix(1+num_frames:end,ii) = visMatrix(1+num_frames:end,ii)-...
        mean(visMatrix(1+num_frames:1+num_frames+6,ii));
end

%%

visAve = mean(visMatrix(winToUse(1):winToUse(2),:));
[val,indices] = sort(visAve,'descend');

% visMatrixAdj = visMatrix(:,[indices(1:50) indices(end-50:end)]);
% holoMatrixAdj = holoMatrix(:,[indices(1:50) indices(end-50:end)]);
%% Perform PCA on the visual data and project the holographic stimulations

% Note: T*W' = training_matrix - mean(training_matrix)
% So, T = training_matrix*W, since W'*W = I
% W can be used to predict the trajectory in PC space of our stim trials
[W_vis, T_vis, eigenvalues_vis] = pca(visMatrix);
holoProjection = (holoMatrix-mean(holoMatrix))*W_vis;

%% Plot the results 
% Plot of percent variance explained
figure(1); clf; 
subplot(1,2,1); hold on;
plot(cumsum(eigenvalues_vis/sum(eigenvalues_vis)),'linewidth',1.5)
set(gca,'fontsize',16)
ylabel('Percenter variance explained')
xlabel('Number of Components')
title('Vis PCA')

figure(2); clf; 
subplot(2,2,1); hold on;
imagesc(visMatrix(1:num_frames,indices)')
caxis([-1 3])
xlim([1 size(visMatrix,1)/2])
ylim([1 size(visMatrix,2)])
title(sprintf('%0.0f degree grating',oriToUse(1)))
set(gca,'fontsize',16)
ylabel('Neuron number')
xlabel('Time frames')
colorbar()

subplot(2,2,2); hold on;
imagesc(visMatrix(1+num_frames:end,indices)')
caxis([-1 3])
xlim([1 size(visMatrix,1)/2])
ylim([1 size(visMatrix,2)])
title(sprintf('%0.0f degree grating',oriToUse(2)))
set(gca,'fontsize',16)
yticks([])
xlabel('Time frames')
colorbar()

subplot(2,2,3); hold on;
imagesc(holoMatrix(1:num_frames,indices)')
caxis([-1 3])
xlim([1 size(visMatrix,1)/2])
ylim([1 size(visMatrix,2)])
title(sprintf('%0.0f holo stim',oriToUse(1)))
set(gca,'fontsize',16)
ylabel('Neuron number')
xlabel('Time frames')
colorbar()

subplot(2,2,4); hold on;
imagesc(holoMatrix(1+num_frames:end,indices)')
caxis([-1 3])
xlim([1 size(visMatrix,1)/2])
ylim([1 size(visMatrix,2)])
title(sprintf('%0.0f holo stim',oriToUse(2)))
set(gca,'fontsize',16)
yticks([])
xlabel('Time frames')
colorbar()

%% Plot the PCs for each stim angle
colorscheme = [0    0.4470    0.7410; 0.8500    0.3250    0.0980];
figure(3); clf; 
for kk = 1:2
    for jj = 1:2
        subplot(2,3,jj+(kk-1)*3); hold on;
        for ii = 1:2
            if jj == 1
                plot(T_vis(1+num_frames*(ii-1):num_frames*ii,kk),'linewidth',1.5,'color',colorscheme(ii,:))
                if kk == 1
                    title('Visual Stim')
                end
            else
                plot(holoProjection(1+num_frames*(ii-1):num_frames*ii,kk),'--','linewidth',1.5,'color',colorscheme(ii,:))
                if kk == 1
                    title('Holo Stim')
                end
            end
        end
        maxVal = max(max(max(T_vis(:,kk))),max(max(holoProjection(:,kk))));
        minVal = min(min(min(T_vis(:,kk))),min(min(holoProjection(:,kk))));
        plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
        plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
        set(gca,'fontsize',16)
        xlabel('Time frames')
        if kk == 1
            ylabel('PC1')
        else
            ylabel('PC2')
        end
        ylim([minVal maxVal])
    end
end

subplot(1,3,3); hold on;
for ii = 1:2
    h_leg(ii) = plot(T_vis(1+num_frames*(ii-1):num_frames*ii,1),T_vis(1+num_frames*(ii-1):num_frames*ii,2),...
        'linewidth',1.5,'color',colorscheme(ii,:));
    plot(T_vis(1+num_frames*(ii-1),1),T_vis(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
    g_leg(ii) = plot(holoProjection(1+num_frames*(ii-1):num_frames*ii,1),...
        holoProjection(1+num_frames*(ii-1):num_frames*ii,2),'--','linewidth',1.5,...
        'color',colorscheme(ii,:));
    
    plot(holoProjection(1+num_frames*(ii-1),1),holoProjection(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
end

set(gca,'fontsize',16)
xlabel('PC1')
ylabel('PC2')

legend([h_leg g_leg],{'180\circ vis','270\circ vis','180\circ holo','270\circ holo'})



%% Perform PCA on the holo data and project the visual stimulations
[W_holo, T_holo, eigenvalues_holo] = pca(holoMatrix);
visProjection = (visMatrix-mean(visMatrix))*W_holo;

figure(1)
subplot(1,2,2); hold on;
plot(cumsum(eigenvalues_holo/sum(eigenvalues_holo)),'linewidth',1.5)
set(gca,'fontsize',16)
ylabel('Percenter variance explained')
xlabel('Number of Components')
title('Holo PCA')

%%
figure(4); clf; 
for kk = 1:2
    for jj = 1:2
        subplot(2,3,jj+(kk-1)*3); hold on;
        for ii = 1:2
            if jj == 1
                plot(T_holo(1+num_frames*(ii-1):num_frames*ii,kk),'--','linewidth',1.5,'color',colorscheme(ii,:))
                if kk == 1
                    title('Holo Stim')
                end
            else
                plot(visProjection(1+num_frames*(ii-1):num_frames*ii,kk),'linewidth',1.5,'color',colorscheme(ii,:))
                if kk == 1
                    title('Vis Stim')
                end
            end
        end
        maxVal = max(max(max(T_holo(:,kk))),max(max(visProjection(:,kk))));
        minVal = min(min(min(T_holo(:,kk))),min(min(visProjection(:,kk))));
        plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
        plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
        set(gca,'fontsize',16)
        xlabel('Time frames')
        if kk == 1
            ylabel('PC1')
        else
            ylabel('PC2')
        end
        ylim([minVal maxVal])
    end
end

subplot(1,3,3); hold on;
for ii = 1:2
    h_leg(ii) = plot(T_holo(1+num_frames*(ii-1):num_frames*ii,1),T_holo(1+num_frames*(ii-1):num_frames*ii,2),...
        '--','linewidth',1.5,'color',colorscheme(ii,:));
    plot(T_holo(1+num_frames*(ii-1),1),T_holo(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
    g_leg(ii) = plot(visProjection(1+num_frames*(ii-1):num_frames*ii,1),...
        visProjection(1+num_frames*(ii-1):num_frames*ii,2),'linewidth',1.5,...
        'color',colorscheme(ii,:));
    
    plot(visProjection(1+num_frames*(ii-1),1),visProjection(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
end

set(gca,'fontsize',16)
xlabel('PC1')
ylabel('PC2')

legend([h_leg g_leg],{'180\circ holo','270\circ holo','180\circ vis','270\circ vis'})

