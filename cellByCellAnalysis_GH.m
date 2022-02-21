%%
% Attempts to recreate the plots using a cell-by-cell analysis, where the
% final average occurs across all cells, as opposed to across ensembles
%
% Unclear how to index into All.out.anal.minDistbyHolo (see line ~96 below)
%
% Also unclear how to baseline All.out.exp.dfData (currently using the
% first six seconds average across relevant trials for the cell; line ~88)
%%
clear; close all; clc;

%%
% Adds all subfolders to the path
restoredefaultpath;
folder = fileparts(which('cellByCellAnalysis_GH.m')); 
addpath(genpath(folder));
rmpath(folder);

%% Loop over each experiment

loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/'; % where ever your files are
loadList_all = oriLoadList_GH('long');

cellData = [];
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
    
    
    %% Only consider valid ensembles from this experiment
    ensIDs = outVars.ensHNumber(outVars.ensemblesToUse);
    
    %% Loop over all ensembles and store the data
    for ii = 1:length(ensIDs)
        
        ensID = ensIDs(ii);
        s = ensID;
        us = unique(All.out.exp.stimID);
        u = us(s);
        
        % Trials to use
        trialsToUse = ismember(All.out.exp.stimID, u) &...
            All.out.exp.lowMotionTrials &...
            All.out.exp.lowRunTrials &...
            All.out.exp.stimSuccessTrial &...
            ismember(All.out.exp.visID, [0 1]);
        
        %% Cells to use
        holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
        
        % tg = All.out.exp.rois{holo}; % this is cells in ensemble
        % fixed finding the cells in the ensembles? 
        % Not used for this analysis
        tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));
        
        offTargetRisks = All.out.anal.offTargetRisk(holo,:);
        ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
        visCells = All.out.anal.pVisR < 0.05; % visually responsive
        
        cellsToUse = ~offTargetRisks & ~ROIinArtifact' & visCells;
        
        %%
        opts.recWinRange = [0.5 1.5];
        recWinSec=opts.recWinRange+All.out.vis.visStart;
        winToUse = round(recWinSec*All.out.info.FR);
        
        tracesHolozdfData = All.out.exp.zdfData(cellsToUse, :, trialsToUse);
        tracesHoloAll = All.out.exp.allData(cellsToUse, :, trialsToUse);
        
        zdfCellResp = mean((mean(tracesHolozdfData(:,winToUse(1):winToUse(end),:),2)),3);
        allCellResp = mean((mean(tracesHoloAll(:,winToUse(1):winToUse(end),:),2)),3);
        
        tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
        dfCellResp = mean((mean(tracesHolodfData(:,winToUse(1):winToUse(end),:),2)),3);
        
        %% Baseline the df data
        baselineEst = mean((mean(tracesHolodfData(:,1:6,:),2)),3);
        dffCellResp = (dfCellResp-baselineEst)./baselineEst;
        
        %% Unclear how to index into All.out.anal.minDistbyHolo
        
        % All.out.exp.ensMaxD has 20 entires but All.out.anal.minDistbyHolo
        % has 21 entries
        % All.out.anal.minDistbyHolo(1,:) is all 0's so using holo seems
        % wrong. So instead, use ensID (or s) from above
        cellDist = All.out.anal.minDistbyHolo(ensID,cellsToUse)';
        
        %% Shared across all cells 
        cellEnsMaxD = All.out.exp.ensMaxD(holo)*ones(sum(cellsToUse),1);
        cellMeanEnsOSI = All.out.exp.meanEnsOSI(holo)*ones(sum(cellsToUse),1);
        cellEnsOSI = All.out.exp.ensOSI(holo)*ones(sum(cellsToUse),1);
        %%
        cellData = [cellData; zdfCellResp allCellResp dfCellResp dffCellResp...
            cellDist cellEnsMaxD cellMeanEnsOSI cellEnsOSI];
        
    end
    
    fprintf('Exp %d out of %d done loading\n',outer_loop, length(loadList_all));
end


%% Plot the data
distBins = [15:25:500];

cellDistDataAve = zeros(1,length(distBins)-1);
cellDistDataMedian = zeros(1,length(distBins)-1);
cellDistDataErr = zeros(1,length(distBins)-1);

ylabels = {'zdf','All','df','dff'};

figure(1); clf;
for ii = 1:4
    
    cellDistData=[];
    for ll = 1:length(distBins)-1
        cellSelector = cellData(:,5)>distBins(ll) & cellData(:,5)<distBins(ll+1);
        
        cellDistDataAve(ll) = nanmean(cellData(cellSelector,ii));
        cellDistDataMedian(ll) = nanmedian(cellData(cellSelector,ii));
        cellDistDataErr(ll) = nanstd(cellData(cellSelector,ii))/sqrt(sum(cellSelector));
    end
    
    subplot(4,2,1+(ii-1)*2)
    plot(cellData(:,5),cellData(:,ii),'.')
    set(gca,'fontsize',16)
    ylabel(ylabels{ii})
    xlabel('Min Dist')
    xlim([0 250])
    
    subplot(4,2,2+(ii-1)*2); hold on;
    leg(1) = plot(distBins(1:end-1),cellDistDataAve,'k.','markersize',16);
    errorbar(distBins(1:end-1),cellDistDataAve,cellDistDataErr,'k','linewidth',1.5)
    leg(2) = plot(distBins(1:end-1),cellDistDataMedian,'.','markersize',16);
    plot([0 distBins(end)],0*[0 distBins(end)],'k--');
    set(gca,'fontsize',16)
    xlabel('Min Dist')
    xlim([0 250])
    
    if ii == 4
        legend(leg,{'Mean','Median'})
    end
end


% Column Information
% cellData = [zdfCellResp allCellResp dfCellResp dffCellResp 
%             cellDist cellEnsMaxD cellMeanEnsOSI cellEnsOSI];

figure(4); clf;

% 1 for zscore; 4 for baseline df
metric_num = 1; 

ensDistThresTight = 400;
ensDistThresLoose = 500;

ensTuningList = {'untuned','mixedTuned','coTuned'};

untunedThres = 0.3;
meanEnsOSIThres = 0.5;
ensOSIThres = 0.5;
ylimMax = 0;
ylimMin = 0;

for jj = 1:length(ensTuningList)
    ensTuning = ensTuningList(jj);
    for ii = 1:2
        
        cellDistData=[];
        for ll = 1:length(distBins)-1
            if strcmp(ensTuning,'untuned')
                if ii == 2
                    cellSelector = cellData(:,5)>distBins(ll) & cellData(:,5)<distBins(ll+1) ...
                        & cellData(:,6) < ensDistThresTight  ...
                        & cellData(:,7) < untunedThres & cellData(:,8) < untunedThres;
                else
                    cellSelector = cellData(:,5)>distBins(ll) & cellData(:,5)<distBins(ll+1) ...
                        & cellData(:,6) > ensDistThresLoose  ...
                        & cellData(:,7) < untunedThres & cellData(:,8) < untunedThres;
                end
            elseif strcmp(ensTuning,'mixedTuned')
                if ii == 2
                    cellSelector = cellData(:,5)>distBins(ll) & cellData(:,5)<distBins(ll+1) ...
                        & cellData(:,6) < ensDistThresTight  ...
                        & cellData(:,7) > meanEnsOSIThres & cellData(:,8) < ensOSIThres;
                else
                    cellSelector = cellData(:,5)>distBins(ll) & cellData(:,5)<distBins(ll+1) ...
                        & cellData(:,6) > ensDistThresLoose  ...
                        & cellData(:,7) > meanEnsOSIThres & cellData(:,8) < ensOSIThres;
                end
            elseif strcmp(ensTuning,'coTuned')
                if ii == 2
                    cellSelector = cellData(:,5)>distBins(ll) & cellData(:,5)<distBins(ll+1) ...
                        & cellData(:,6) < ensDistThresTight  ...
                        & cellData(:,7) > meanEnsOSIThres & cellData(:,8) > ensOSIThres;
                else
                    cellSelector = cellData(:,5)>distBins(ll) & cellData(:,5)<distBins(ll+1) ...
                        & cellData(:,6) > ensDistThresLoose  ...
                        & cellData(:,7) > meanEnsOSIThres & cellData(:,8) > ensOSIThres;
                end
            end
            
            cellDistDataAve(ll) = nanmean(cellData(cellSelector,metric_num));
            cellDistDataMedian(ll) = nanmedian(cellData(cellSelector,metric_num));
            cellDistDataErr(ll) = nanstd(cellData(cellSelector,metric_num))/sqrt(sum(cellSelector));
        end
        
        subplot(2,3,ii+(ii-1)*2+jj-1); hold on;
        leg(1) = plot(distBins(1:end-1),cellDistDataAve,'k.','markersize',16);
        errorbar(distBins(1:end-1),cellDistDataAve,cellDistDataErr,'k','linewidth',1.5)
        leg(2) = plot(distBins(1:end-1),cellDistDataMedian,'.','markersize',16);
        set(gca,'fontsize',16)
        xlabel('Min Dist')
        xlim([0 150])
        plot([0 distBins(end)],0*[0 distBins(end)],'k--');
        
        ylimMax = max(max(ylim),ylimMax);
        ylimMin = min(min(ylim),ylimMin);
    end 
end


subplot(2,3,1)
title(sprintf('Untuned \n ensOSI < %.1f, meanEnsOSI < %.1f',untunedThres,untunedThres))
if metric_num==1 
    ylabel(sprintf('Max Ens Dist > 500 \n z-scored df'))
elseif metric_num==4
    ylabel(sprintf('Max Ens Dist > 500 \n (f-f0)/f0'))
end
subplot(2,3,2)
title(sprintf('Mixed Tuned \n meanEnsOSI > %.1f, ensOSI < %.1f',meanEnsOSIThres, ensOSIThres))
subplot(2,3,3)
title(sprintf('Cotuned \n meanEnsOSI > %.1f, ensOSI > %.1f',meanEnsOSIThres,ensOSIThres))
subplot(2,3,4)
if metric_num==1 
    ylabel(sprintf('Max Ens Dist < 400 \n z-scored df'))
elseif metric_num==4
    ylabel(sprintf('Max Ens Dist < 400 \n (f-f0)/f0'))
end


for ii = 1:6
    subplot(2,3,ii)
   ylim([ylimMin ylimMax]) 
end
