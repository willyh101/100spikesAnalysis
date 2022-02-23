clear; close all; clc;

%% Load data
% Options: default, all, stimmedAndOT, stimmed
includedCells = 'default';

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

ensIDs  = zeros(1,length(oriToUse));
ensOri = zeros(1,length(oriToUse));
for ii = 1:length(oriToUse)
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
% orisToUseVis = [nan oriToUse(2)]; % copy from exp  section

if isnan(orisToUseVis(1))
    trialsToUse = (ismember(oris, orisToUseVis)+isnan(oris))>0;
else
    trialsToUse = ismember(oris, orisToUseVis);
end

% cells to use is copied from above
% result is cells x time x trials
tracesVis = visData(cellsToUse, :, trialsToUse);
trialAngle = oris(trialsToUse);

%%

%% Split the trials into two (training/testing), and average across them
visTestMatrix = zeros(num_frames*2,size(tracesVis,1));
visMatrix = zeros(num_frames*2,size(tracesVis,1));
% average across trials
for ii = 1:length(ensIDs)
       
    if isnan(orisToUseVis(ii))
        tempArray = find(isnan(trialAngle)==1);
    else
        tempArray = find(trialAngle==oriToUse(ii));
    end
    totalTrials = length(tempArray);
    if totalTrials>10
        firstHalf = floor(totalTrials/2);
        tempArray = tempArray(randperm(totalTrials));
        
        trialAveVis = mean(tracesVis(:,:,tempArray(1:firstHalf)),3);
        visMatrix(1+(ii-1)*num_frames:ii*num_frames,:) = trialAveVis';
        
        trialAveVis = mean(tracesVis(:,:,tempArray(1+firstHalf:end)),3);
        visTestMatrix(1+(ii-1)*num_frames:ii*num_frames,:) = trialAveVis';
    else
        error('Probably not enough trials to proceed')
    end
    
end

%% Optional pre-processing (zero-mean before stimulus across trials) 
for ii = 1:size(visMatrix,2)
    for jj = 1:length(ensIDs)
        visMatrix(1+num_frames*(jj-1):num_frames*jj,ii) ...
            = visMatrix(1+num_frames*(jj-1):num_frames*jj,ii)-...
            mean(visMatrix(1+num_frames*(jj-1):num_frames*(jj-1)+6,ii));
        
        visTestMatrix(1+num_frames*(jj-1):num_frames*jj,ii) ...
            = visTestMatrix(1+num_frames*(jj-1):num_frames*jj,ii)-...
            mean(visTestMatrix(1+num_frames*(jj-1):num_frames*(jj-1)+6,ii));
    end
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
visProjection = (visTestMatrix-mean(visTestMatrix))*W_vis;

%%
colorscheme = [0    0.4470    0.7410; 0.8500    0.3250    0.0980];
%%

figure(4); clf; 
subplot(1,3,1); hold on;
for ii = 1:2
    plot(T_vis(1+num_frames*(ii-1):num_frames*ii,1),'-','linewidth',1.5,'color',colorscheme(ii,:))
    plot(visProjection(1+num_frames*(ii-1):num_frames*ii,1),'--','linewidth',1.5,'color',colorscheme(ii,:))
end
maxVal = max(max(max(T_vis(:,1))),max(max(visProjection(:,1))));
minVal = min(min(min(T_vis(:,1))),min(min(visProjection(:,1))));
plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
set(gca,'fontsize',16)
xlabel('Time frames')
ylabel('PC1')
ylim([minVal maxVal])

subplot(1,3,2); hold on;
for ii = 1:2
    plot(T_vis(1+num_frames*(ii-1):num_frames*ii,2),'-','linewidth',1.5,...
        'color',colorscheme(ii,:))
    
    plot(visProjection(1+num_frames*(ii-1):num_frames*ii,2),'--','linewidth',1.5,...
        'color',colorscheme(ii,:))
end
maxVal = max(max(max(T_vis(:,2))),max(max(visProjection(:,2))));
minVal = min(min(min(T_vis(:,2))),min(min(visProjection(:,2))));
set(gca,'fontsize',16)
xlabel('Time')
ylabel('PC2')
plot(winToUse(1)+0*[minVal maxVal],[minVal maxVal],'k--')
plot(winToUse(2)+0*[minVal maxVal],[minVal maxVal],'k--')
ylim([minVal maxVal])

subplot(1,3,3); hold on;
for ii = 1:2
    g_leg(ii) = plot(T_vis(1+num_frames*(ii-1):num_frames*ii,1),T_vis(1+num_frames*(ii-1):num_frames*ii,2),...
        '-','linewidth',1.5,'color',colorscheme(ii,:));
    plot(T_vis(1+num_frames*(ii-1),1),T_vis(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
    
    h_leg(ii) = plot(visProjection(1+num_frames*(ii-1):num_frames*ii,1),visProjection(1+num_frames*(ii-1):num_frames*ii,2),...
        '--','linewidth',1.5,'color',colorscheme(ii,:));
    plot(visProjection(1+num_frames*(ii-1),1),visProjection(1+num_frames*(ii-1),2),'k.','markersize',16,'linewidth',1.5)
    
end

set(gca,'fontsize',16)
xlabel('PC1')
ylabel('PC2')

legend([g_leg h_leg],{strcat(sprintf('%0.0f',oriToUse(1)),'\circ vis'),...
    strcat(sprintf('%0.0f',oriToUse(2)),'\circ vis'),...
    strcat(sprintf('%0.0f',oriToUse(1)),'\circ vis'),...
    strcat(sprintf('%0.0f',oriToUse(2)),'\circ vis')})

%%

figure(1); clf; hold on;
plot(cumsum(eigenvalues_vis/sum(eigenvalues_vis)),'linewidth',1.5)
set(gca,'fontsize',16)
ylabel('Percenter variance explained')
xlabel('Number of Components')
title('Vis PCA')

participation_ratio = sum(eigenvalues_vis)/sum(eigenvalues_vis.^2)


%%
figure(2); clf; 
for jj = 1:2
    for ii = 1:2
        subplot(2,2,ii+(jj-1)*2); hold on;
        if jj == 1
            imagesc(visMatrix(1+(ii-1)*num_frames:ii*num_frames,indices)')
        else
            imagesc(visTestMatrix(1+(ii-1)*num_frames:ii*num_frames,indices)')
        end
        caxis([0 1])
        xlim([1 size(visMatrix,1)/2])
        ylim([1 size(visMatrix,2)])
        title(sprintf('%0.0f degree grating',oriToUse(ii)))
        set(gca,'fontsize',16)
        ylabel('Neuron number')
        xlabel('Time frames')
        colorbar()
    end
end
%%

figure(11); clf; hold on;
aveResponse = mean(visMatrix,2);
aveTestResponse = mean(visTestMatrix,2);
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
