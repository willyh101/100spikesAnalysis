%%
% Recreates the key plots of the paper (averaged over ensembles), and 
% reproducing them using a cell-by-cell analysis 
%%
clear; close all; clc;


%%
% Adds all subfolders to the path
restoredefaultpath;
folder = fileparts(which('cellByCellAnalysis_GH.m')); 
addpath(genpath(folder));
rmpath(folder);

%% Specify data location and loadList

% loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/'; % where ever your files are
loadPath = '/Users/greghandy/Research_Local_v2/'; % where ever your files are

% Options: short (for debugging), all (all avaliable outfiles), used (data
% files currently being used here; created for speed)
loadList_all = oriLoadList_GH('used');

%% Loop over each experiment
ensNum = 1;
expNum = 1;
cellData = [];
expData=[];
cellID=0;
%%
for outer_loop = 1:length(loadList_all) %10 for 60 offTarget
    
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
    %% Additional pre-processing: 
    % Remove cells from trials where they were offTargets/targets within some
    % buffer zone
    % Key parameters: length of buffer zone and whether to remove targets
    % or offTargets
    % If the total trial count for these cells drops below 10, they are
    % removed entirely
    
    % length(All.out.anal.targets) % total unique targets
    unID = unique(All.out.exp.stimID);
    IDshift = unID(2)-1;
    totalCells = size(All.out.exp.dfData,1);
    totalTrials = size(All.out.exp.stimID,2);
    
    % Everyone starts being avaliable for analysis
    bufferZone = 1;
    lastTimeStimmed = bufferZone+zeros(totalCells,1);
    % Loop through all trials sequentially
    for ii = 1:totalTrials
        
        % Remove these cells from this trial
        doNotInclude = lastTimeStimmed<bufferZone;
        All.out.exp.dfData(doNotInclude,:,ii) = nan;
        
        % Update the lastTimeStimmed list
        % Figure out who is currently stimmed
        holoIndex = All.out.exp.stimID(ii)-IDshift;
        if holoIndex > 0
            % If you only want to remove targets, uncomment this line
            currTargets=find(All.out.anal.offTargetRisk(holoIndex,:)==1);
        else
            currTargets=[];
        end
        lastTimeStimmed = lastTimeStimmed + 1;
        % Reset the count of those just stimmed
        lastTimeStimmed(currTargets) = 0;
    end
    
    %% Only consider valid ensembles from this experiment
    ensIDs = outVars.ensHNumber(outVars.ensemblesToUse);
    cellsToUse = [];
    tempTrial = 0;
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
            All.out.exp.visID == 1;                
        tempTrial = tempTrial + sum(trialsToUse);
        %% Cells to use
        ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
        cellsToUse = ~ROIinArtifact' & ~All.out.anal.cellsToExclude;
        
        %% pVis value
        pVis = All.out.anal.pVisR;
        
        %% OffTargets and targets
        holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
        % These are the cell indices for those in ensemble
        tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));
        tgCell = zeros(length(cellsToUse),1);
        tgCell(tg_adj) = true;
        
        offTargetRisks = All.out.anal.offTargetRisk(holo,:);
        %% Baseline the df data
        % Options: dfData, zdfData, allData
        recWinSec=opts.recWinRange+All.out.exp.visStart;
        winToUse = round(recWinSec*All.out.info.FR);
        bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
        tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
        dfCellResp = nanmean((nanmean(tracesHolodfData(:,winToUse(1):winToUse(end),:),2)),3);
        baselineEst = nanmean((nanmean(tracesHolodfData(:,bwinToUse(1):bwinToUse(2),:),2)),3);   
        dffCellResp = (dfCellResp-baselineEst);
        
        % Remove cells that don't have enough (10) trials
        tempTotalTrials = sum(trialsToUse);
        badCells = (tempTotalTrials-sum(isnan(tracesHolodfData(:,1,:)),3))<10;
        dffCellResp(badCells) = nan;
                
        %% Control response (using the no stim epochs)
        % Method 1: 
        stimIDs = unique(All.out.exp.stimID);        
        trialsToUse = ismember(All.out.exp.stimID, stimIDs(1)) &...
            All.out.exp.lowMotionTrials &...
            All.out.exp.lowRunTrials & All.out.exp.visID==1;
        tempTotalTrials = sum(trialsToUse);
        
        recWinSec=opts.recWinRange+All.out.exp.visStart;
        winToUse = round(recWinSec*All.out.info.FR);
        bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
        tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);        
        dfCellResp = nanmean((nanmean(tracesHolodfData(:,winToUse(1):winToUse(end),:),2)),3);
        baselineEst = nanmean((nanmean(tracesHolodfData(:,bwinToUse(1):bwinToUse(2),:),2)),3);
        dffCellCtlResp_noStim = (dfCellResp-baselineEst);
        
        % Remove cells that don't have enough (10) trials
        badCells = (tempTotalTrials-sum(isnan(tracesHolodfData(:,1,:)),3))<10;
        dffCellCtlResp_noStim(badCells) = nan;
            
        if unique(All.out.anal.minDistbyHolo(1,:)) ~=0
            error('Something wrong with minDistbyHolo');
        end
        
        %% Careful indexing into All.out.anal.minDistbyHolo
        cellDist = All.out.anal.minDistbyHolo(ensID,cellsToUse)';
        
        %% These are not ensemble-dependent
        cellPO = All.out.anal.prefOri(cellsToUse)';
        cellOSI = All.out.anal.osi(cellsToUse)';
        
        %% Shared across all cells 
        cellEnsMaxD = All.out.exp.ensMaxD(holo)*ones(sum(cellsToUse),1);
        cellMeanEnsOSI = All.out.exp.meanEnsOSI(holo)*ones(sum(cellsToUse),1);
        cellEnsOSI = All.out.exp.ensOSI(holo)*ones(sum(cellsToUse),1);
        cellEnsPO = All.out.exp.ensPO(holo)*ones(sum(cellsToUse),1);
                
        %% Store the data
        if ii==1
            IDStart = cellID(end)+1;
            cellID = (IDStart:(IDStart-1)+sum(cellsToUse))';
        end
               
        cellData = [cellData; dffCellResp cellDist cellPO cellOSI ...
            cellEnsMaxD cellMeanEnsOSI cellEnsOSI pVis(cellsToUse)' offTargetRisks(cellsToUse)'...
            tgCell(cellsToUse) ensNum*ones(sum(cellsToUse),1) expNum*ones(sum(cellsToUse),1) cellEnsPO,...
            dffCellCtlResp_noStim cellID];
        ensNum = ensNum + 1;
    end
    
    if ~isempty(cellsToUse)
        expNum = expNum +1;
    end
    
    %%
    
    fprintf('Exp %d out of %d done loading\n',outer_loop, length(loadList_all));
end

%% Turn matrix into a table (easier to access specific columns by name)

cellTable = array2table(cellData,'VariableNames',{'dff','cellDist','cellPO','cellOSI','cellEnsMaxD',...
    'cellMeanEnsOSI','cellEnsOSI','visP','offTarget','tgCell','ensNum','expNum','ensPO',...
    'ctlResp_noStim','cellID'});

%% Iso vs. ortho calculation
oriVals = [NaN 0:45:315];
cellOris = oriVals(cellTable.cellPO)';
cellOrisDiff = abs(cellOris-cellTable.ensPO);
cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
cellOrisDiff(cellOrisDiff==135)=45;
cellOrisDiff(cellOrisDiff==180)=0;
cellTable.cellOrisDiff = cellOrisDiff;

%% Cell conditions used in the functions
cellCond = cellTable.offTarget==0;
cellCondTuned = cellTable.offTarget==0 & cellTable.visP<0.05 & cellTable.cellOSI > 0.25;

%% Figure 3: min distance plot
Fig3(cellTable,cellCond)
Fig3_cbc(cellTable,cellCond)
Fig3_noStimCtl(cellTable)

%% Percent of cells suppressed vs. activated across all ensembles
allResp = cellTable.dff(cellCond);
allResp(isnan(allResp))=[];
temp1=sum(allResp<0)/sum(~isnan(allResp));
temp2=sum(allResp>0)/sum(~isnan(allResp));

figure(7773); clf; hold on
bar([temp1;temp2]*100)
set(gca,'fontsize',16)
xticks([1 2])
xticklabels({'Suppressed','Activated'})
ylabel('Percent of cells')
plot([0 3],50+[0 0],'k--','linewidth',1.5)
xlim([0 3])

%% Figure Space: Effect of the spread of ensemble
ensResp = FigSpace(cellTable,cellCond);

%% Figure 5: Iso vs. ortho
Fig5(cellTable,cellCondTuned);
Fig5_cbc(cellTable,cellCondTuned);

%% Figure 6: Tight co-tuned investigation
Fig6(cellTable,cellCond);
Fig6_cbc(cellTable,cellCond);
Fig6IsoOrtho(cellTable,cellCondTuned)

%% Figure plotting the percent activated/suppressed as a function of dist
% Third input is threshold value
FigPercentAct(cellTable,cellCond,0)

%% Target cell statistics
FigTargets(cellTable)

%% Figure illustrating the total number of cells suppressed
clc;
FigCellByCell(cellTable,cellCond)

%%

totalExps = cellTable.expNum(end);
grayCtl_noStim = [];
for ii = 1:totalExps
    tempIndices = find(cellTable.expNum == ii);
    tempEns = unique(cellTable.ensNum(tempIndices));
    % We might be getting weird edge effects for larger distances
    cellSelector = cellTable.ensNum == tempEns(1) & cellCond;
    
    grayCtl_noStim(ii) = nanmean(cellTable.ctlResp_noStim(cellSelector));
    perPost_noStim(ii) = sum(cellTable.ctlResp_noStim(cellSelector)<0)/sum(cellSelector);
end

ensAve_noStim = nanmean(grayCtl_noStim);
ensSTD_noStim = nanstd(grayCtl_noStim)/sqrt(length(grayCtl_noStim));

figure(324); clf;
subplot(1,2,1); hold on
plot(1+0.1*rand(length(grayCtl_noStim),1),grayCtl_noStim,'.','markersize',12)
errorbar(1.05,ensAve_noStim,ensSTD_noStim,'linewidth',2,'color','k')
plot([0.5 1.5],0*[0.5 1.5],'k--','markersize',16)
set(gca,'fontsize',16)
xticks([1])
xticklabels({'Gray screen (no stim)'})
ylabel('dF/F')

subplot(1,2,2); hold on;
boxplot(perPost_noStim*100)
% xlim([0.5 2.5])
xticks([1])
xticklabels({'Gray (no stim)'})
set(gca,'fontsize',16)
ylabel('% of cells suppressed')
