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
        %(All.out.exp.visID == 1 |  All.out.exp.visID == 0); % not sure what this is (1 is Gray; 0 is nothing?) 
        %ismember(All.out.exp.visID, [0 1]);
                
        tempTrial = tempTrial + sum(trialsToUse);
        %% Cells to use
        holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
        
        % tg = All.out.exp.rois{holo}; % these ROIs for cells in ensemble
        % These are the cell indices for those in ensemble
        tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));
        
        offTargetRisks = All.out.anal.offTargetRisk(holo,:);
        offTargetRisksOnePlane = All.out.anal.offTargetRiskOnePlane(holo,:);
        ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
        pVis = All.out.anal.pVisR;
        % visCells = All.out.anal.pVisR < 0.05; % visually responsive
        % cellsToUse = ~offTargetRisks & ~ROIinArtifact' & visCells;
        cellsToUse = ~ROIinArtifact' & ~All.out.anal.cellsToExclude;
                    
        tgCell = zeros(length(cellsToUse),1);
        tgCell(tg_adj) = true;
        numOfTgs(ensNum) = length(tg_adj);
        
        radialDist = All.out.anal.radialDistToStim(holo,cellsToUse)';
        %% Baseline the df data
        % Options: dfData, zdfData, allData
        
        % Method 1: 
        recWinSec=opts.recWinRange+All.out.exp.visStart;
        winToUse = round(recWinSec*All.out.info.FR);
        bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
        tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
        dfCellResp = mean((mean(tracesHolodfData(:,winToUse(1):winToUse(end),:),2)),3);
        baselineEst = mean((mean(tracesHolodfData(:,bwinToUse(1):bwinToUse(2),:),2)),3);
        dffCellResp2 = (dfCellResp-baselineEst);
         
        % Method 2: use results from anal 
        tempResp = squeeze(All.out.anal.respMat(s,cellsToUse));
        tempBase = squeeze(All.out.anal.baseMat(s,cellsToUse));
        dffCellResp = (tempResp-tempBase)';
        
        % Methods should be identical
        if norm(dffCellResp-dffCellResp2)~=0
           plot(dffCellResp,dffCellResp2)
           'here';
        end

        %% Control response (using the gray screen epochs from vis sections)
        recWinSec=opts.recWinRange+All.out.vis.visStart;
        winToUse = round(recWinSec*All.out.info.FR);
        bwinToUse = max(floor([0 All.out.vis.visStart]*All.out.info.FR),[1 1]);
        
        % visID == 1 should be a gray screen
        visTrialsToUse = All.out.vis.visID == 1; 
        % Including these restrictions results in no usable trials for experiment 3
        % & All.out.vis.lowRunTrials & All.out.vis.lowMotionTrials;
             
        tempCtlResp = nanmean(All.out.vis.rdData(:,visTrialsToUse),2);
        tempCtlBase = nanmean(All.out.vis.bdata(:,visTrialsToUse),2);
        dffCellCtlResp_vis = (tempCtlResp-tempCtlBase);
        
        %% Control response (using the no stim epochs)
%         stimIDs = unique(All.out.exp.stimID);

        % Method 1: 
        stimIDs = unique(All.out.exp.stimID);
        trialsToUse = ismember(All.out.exp.stimID, stimIDs(1)) &...
            All.out.exp.lowMotionTrials &...
            All.out.exp.lowRunTrials & All.out.exp.visID==1; 
        
        recWinSec=opts.recWinRange+All.out.exp.visStart;
        winToUse = round(recWinSec*All.out.info.FR);
        bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
        tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
        dfCellResp = mean((mean(tracesHolodfData(:,winToUse(1):winToUse(end),:),2)),3);
        baselineEst = mean((mean(tracesHolodfData(:,bwinToUse(1):bwinToUse(2),:),2)),3);
        dffCellCtlResp_v3 = (dfCellResp-baselineEst);

        % Method 2:
        tempCtlResp = squeeze(All.out.anal.respMat(1,cellsToUse));
        tempCtlBase = squeeze(All.out.anal.baseMat(1,cellsToUse));
        dffCellCtlResp_noStim = (tempCtlResp-tempCtlBase)';
        
        %%
%         tgCell_adj = tgCell(cellsToUse);
%         temp = tracesHolodfData(tgCell_adj==1,:,:);
%         tempAve = mean(temp,3)-baselineEst(tgCell_adj==1);
%         tempAve2 = mean(tempAve',2);
%         
%         tAxis = [0:1:size(tempAve2)-1]/All.out.info.FR;
%         
%         figure(1); clf; hold on;
%         plot(tAxis,tempAve','color',[0 0.447 0.741 0.25])
%         plot(tAxis,tempAve2,'color',[0 0.447 0.741],'linewidth',2)
%         plot([0 24]/All.out.info.FR, 0*[0 24],'k--')
%         title('Target cells (no stim trial ave)')
%         xlabel('Time (seconds)')
%         ylabel('\DeltaF/F')
%         set(gca,'fontsize',16)
%         xlim([0 24]/All.out.info.FR)
%         'here';
        
        %%
        
        if norm(dffCellCtlResp_v3-dffCellCtlResp_noStim)~=0
            norm(dffCellCtlResp_v3-dffCellCtlResp_noStim)
            'here'
        end
        
        if unique(All.out.anal.minDistbyHolo(1,:)) ~=0
            'here';
        end
         
        
        %% Careful indexing into All.out.anal.minDistbyHolo
        cellDist = All.out.anal.minDistbyHolo(ensID,cellsToUse)';
        cellPO = All.out.anal.prefOri(cellsToUse)';
        cellOSI = All.out.anal.osi(cellsToUse)';
        
        %% Shared across all cells 
        cellEnsMaxD = All.out.exp.ensMaxD(holo)*ones(sum(cellsToUse),1);
        cellMeanEnsOSI = All.out.exp.meanEnsOSI(holo)*ones(sum(cellsToUse),1);
        cellEnsOSI = All.out.exp.ensOSI(holo)*ones(sum(cellsToUse),1);
        cellEnsPO = All.out.exp.ensPO(holo)*ones(sum(cellsToUse),1);
        %%
        
%         tempIndices = find(dffCellResp>2 & cellDist<20 & ~offTargetRisks(cellsToUse)');
        tempIndices = find(dffCellResp>2 & cellDist>180 & ~offTargetRisks(cellsToUse)');
        if ~isempty(tempIndices) && 0
            tempIndices %(664)
            %%
%             tempIndices = tg_adj();
            tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
            
            recWinSec=opts.recWinRange+All.out.exp.visStart;
            winToUse = round(recWinSec*All.out.info.FR);
            bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
            baselineEst = squeeze(mean(tracesHolodfData(tempIndices(1),bwinToUse(1):bwinToUse(2),:)));
            
            figure(131432499); clf; hold on; 
%             subplot(2,length(tg_adj),1); hold on;
            baselinedResp = zeros(12,24);
            for gg = 1:size(tracesHolodfData,3)
                baselinedResp(gg,:) = tracesHolodfData(tempIndices(1),:,gg)-baselineEst(gg);
                plot(baselinedResp(gg,:),'color',[0 0.447 0.741 0.3])
            end
            
            plot(mean(baselinedResp),'linewidth',2,'color',[0 0.447 0.741])
            plot(round(All.out.exp.visStart*All.out.info.FR)+[0 0],[-1 7],'k--')
            plot(round((All.out.exp.visStart+1)*All.out.info.FR)+[0 0],[-1 7],'k--')
          
            tempResp = mean(baselinedResp);
            mean(tempResp(winToUse(1):winToUse(2)))
            dffCellResp(tempIndices(1))
%             ylim([-1 5])
            title(sprintf('Dist to Stim: %.2f',cellDist(tempIndices(1))))

            %%
        end
        %%
        if  ~isempty(tempIndices) && 0 % expNum>15 &&  0
%             figure(131432499); clf; 
        tracesHolodfData = All.out.exp.dfData(:, :, trialsToUse);
        recWinSec=opts.recWinRange+All.out.exp.visStart;
        winToUse = round(recWinSec*All.out.info.FR);
        bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
        figure(131432499); hold on;
        baselinedResp=[];
        for stimLoop = 1:length(tg_adj)
            
            baselineEst = squeeze(mean(tracesHolodfData(tg_adj(stimLoop),bwinToUse(1):bwinToUse(2),:)));
            subplot(2,length(tg_adj),stimLoop+length(tg_adj)); hold on;
            for gg = 1:size(tracesHolodfData,3)
                baselinedResp(gg,:) = tracesHolodfData(tg_adj(stimLoop),:,gg)-baselineEst(gg);
                plot(baselinedResp(gg,:),'color',[0 0 0 0.3])
            end
            plot(mean(baselinedResp),'k','linewidth',2)
            plot(round(All.out.exp.visStart*All.out.info.FR)+[0 0],[-1 5],'k--')
            plot(round((All.out.exp.visStart+1)*All.out.info.FR)+[0 0],[-1 5],'k--')
            title(sprintf('Dist to Stim: %.2f',All.out.anal.minDistbyHolo(ensID,tg_adj(stimLoop))))
        end
        'here'
        end
        
        %%
        
        cellData = [cellData; dffCellResp cellDist cellPO cellOSI ...
            cellEnsMaxD cellMeanEnsOSI cellEnsOSI pVis(cellsToUse)' offTargetRisks(cellsToUse)'...
            tgCell(cellsToUse) ensNum*ones(sum(cellsToUse),1) expNum*ones(sum(cellsToUse),1) cellEnsPO,...
            dffCellCtlResp_vis(cellsToUse) dffCellCtlResp_noStim offTargetRisksOnePlane(cellsToUse)' ...
            radialDist];
        ensNum = ensNum + 1;
    end
    
    if ~isempty(cellsToUse)
        expNum = expNum +1;
    end
    
    
    %%
    totalEns = length(All.out.exp.holoTargets);
    allTGList = [];
    for hh = 1:totalEns
        allTGList = [allTGList All.out.exp.holoTargets{hh}(~isnan(All.out.exp.holoTargets{hh}))];
    end
    allTGList = unique(allTGList);
    
    anytgCell = zeros(length(cellsToUse),1);
    anytgCell(allTGList) = true;
    
    expData = [expData; (expNum-1)*ones(sum(cellsToUse),1) anytgCell(cellsToUse)];
    %%
    
    
    fprintf('Exp %d out of %d done loading\n',outer_loop, length(loadList_all));
end

%% Turn matrix into a table (easier to access specific columns by name)

cellTable = array2table(cellData,'VariableNames',{'dff','cellDist','cellPO','cellOSI','cellEnsMaxD',...
    'cellMeanEnsOSI','cellEnsOSI','visP','offTarget','tgCell','ensNum','expNum','ensPO','ctlResp_vis',...
    'ctlResp_noStim','offTargetOnePlane','radialDist'});

%% Iso vs. ortho calculation
oriVals = [NaN 0:45:315];
cellOris = oriVals(cellTable.cellPO)';
cellOrisDiff = abs(cellOris-cellTable.ensPO);
cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
cellOrisDiff(cellOrisDiff==135)=45;
cellOrisDiff(cellOrisDiff==180)=0;
cellTable.cellOrisDiff = cellOrisDiff;

%% New column: label cells that were offTarget for any experiment
anyOffTarget = [];
anyTarget = [];
for ee = 1:cellTable.expNum(end)
    expList = unique(cellTable.ensNum(cellTable.expNum==ee));
    numCells = sum(cellTable.ensNum==expList(1));
    anyOffTargetTemp = zeros(numCells,1);
    anyTargetTemp = zeros(numCells,1);
    for ii = 1:length(expList)
        temp1 = cellTable.offTarget(cellTable.ensNum==expList(ii))==1;
%         temp1 = cellTable.radialDist(cellTable.ensNum==expList(ii))<=10 &...
%             cellTable.offTargetOnePlane(cellTable.ensNum==expList(ii))==1;
%         temp3 = cellTable.radialDist(cellTable.ensNum==expList(ii))<=10 &...
%             cellTable.offTargetOnePlane(cellTable.ensNum==expList(ii))==0;
        anyOffTargetTemp = temp1+anyOffTargetTemp;
        
        temp2 = cellTable.tgCell(cellTable.ensNum==expList(ii))==1;
        anyTargetTemp = temp2+anyTargetTemp;
        
    end
    anyOffTargetTemp = sign(anyOffTargetTemp);
    anyOffTarget = [anyOffTarget; repmat(anyOffTargetTemp,length(expList),1)];
    
    anyTargetTemp = sign(anyTargetTemp);  
    anyTarget = [anyTarget; repmat(anyTargetTemp,length(expList),1)];
end
cellTable.anyOffTarget= anyOffTarget;
cellTable.anyTarget= anyTarget;

%%

totalTgs = 0;
tempCol = [];
for ee = 1:cellTable.expNum(end)
    ensList = unique(cellTable.ensNum(cellTable.expNum==ee));
    
    rowSelector = expData(:,1)==ee;
    tempCol = [tempCol; repmat(expData(rowSelector,2),length(ensList),1)];
    
    totalTgs = totalTgs + sum(cellTable.anyTarget(cellTable.ensNum==ensList(1))==1);
end
cellTable.ALLTargets = tempCol;

%% Cell conditions used in the functions
cellCond = cellTable.offTarget==0;
cellCond_v3 = cellTable.anyOffTarget==0;
cellCond_v2 = cellTable.offTarget==0 & cellTable.ALLTargets == 0;
cellCondTuned = cellTable.offTarget==0 & cellTable.visP<0.05 & cellTable.cellOSI > 0.25;
cellCondTuned_v2 = cellTable.offTarget==0 & cellTable.visP<0.05 & cellTable.cellOSI > 0.25 & cellTable.ALLTargets == 0;


%% Figure 3: min distance plot

Fig3(cellTable,cellCond)

%%
% cellCond = cellTable.offTargetOnePlane==1;
Fig3_cbc(cellTable,cellCond)

%%
Fig_noStimCtl(cellTable,expData)
% Fig3_ctl(cellTable)

%% Figure Space: Effect of the spread of ensemble
ensResp = FigSpace(cellTable,cellCond);

%% Figure 5: Iso vs. ortho

Fig5(cellTable,cellCondTuned);
Fig5_cbc(cellTable,cellCondTuned);

%% Figure 6: Tight co-tuned investigation
Fig6(cellTable,cellCond_v2)

%%
Fig6_cbc(cellTable,cellCond_v2);

%%
Fig6IsoOrtho(cellTable,cellCondTuned_v2)

%% Figure plotting the percent activated/suppressed as a function of dist
% Third input is threshold value
FigPercentAct(cellTable,cellCond,0)

%% Target cell statistics
FigTargets(cellTable)

%% Figure illustrating the total number of cells suppressed
clc;
FigCellByCell(cellTable,cellCond)

%%

% ensThreshs high value: neurons share similar POs
ensSelector = cellTable.cellMeanEnsOSI>0.5;
ensOSI = nan(160,1);
ensSpread = nan(160,1);
ensVal = nan(160,1);
for ii = 1:160
    cellSelectorEns = cellTable.ensNum == ii;
    
    cellSelector = cellSelectorEns & cellCond & ensSelector;
    
    if sum(cellSelector)>0
        ensOSI(ii) = unique(cellTable.cellEnsOSI(cellSelector));
        ensSpread(ii) = unique(cellTable.cellEnsMaxD(cellSelector));
        ensVal(ii) = nanmean(cellTable.dff(cellSelector));
    end
end

figure(31); clf; hold on;
scatter(ensSpread,ensOSI,60,ensVal,'filled')


plot(400+[0 0],[0 1],'k--')
plot(500+[0 0],[0 1],'k--')
plot([200 650],0.7 + [0 0],'k--')
plot([200 650],0.3 + [0 0],'k--')
caxis([-0.05 0.03])

colorbar
set(gca,'fontsize',16)
xlabel('Ensemble Spread')
ylabel('Ensemble tuning')

%%

dffBins = [-1:0.1:1];
plotDffBins = dffBins(1:end-1)+0.025;
cellDistAve  = zeros(length(dffBins)-1,1);
cellDistSTDErr= zeros(length(dffBins)-1,1);
numCells= zeros(length(dffBins)-1,1);
for ii = 1:length(dffBins)-1
    cellDFF = cellTable.dff>dffBins(ii) & cellTable.dff<dffBins(ii+1) ;
    cellSelector =  cellCond & cellDFF;


    cellDistAve(ii) = nanmean(cellTable.cellDist(cellSelector));
    cellDistSTDErr(ii) = nanstd(cellTable.cellDist(cellSelector))/sqrt(sum(cellSelector));
    numCells(ii) = sum(cellSelector);
end

figure(132421); clf; hold on;
plot(plotDffBins,cellDistAve,'.','color',[0 0.447 0.741],'markersize',16)

errorbar(plotDffBins,cellDistAve,cellDistSTDErr,'color',[0 0.447 0.741])
set(gca,'fontsize',16)
xlabel('\DeltaF/F')
ylabel('Average min distance to ensemble')
xlim([-1.1 1.1])

%%
actPerPlot(cellTable,cellCond_v2)

%%

totalExps = 18;
grayCtl_noStim = [];
grayCtl_vis = [];
for ii = 1:totalExps
    tempIndices = find(cellTable.expNum == ii);
    tempEns = unique(cellTable.ensNum(tempIndices));
    % We might be getting weird edge effects for larger distances
    cellSelector = cellTable.ensNum == tempEns(1) & cellCond;
    
    grayCtl_noStim(ii) = nanmean(cellTable.ctlResp_noStim(cellSelector));
    grayCtl_vis(ii) = nanmean(cellTable.ctlResp_vis(cellSelector));
    
    perPost_noStim(ii) = sum(cellTable.ctlResp_noStim(cellSelector)>0)/sum(cellSelector);
    perPost_vis(ii) = sum(cellTable.ctlResp_vis(cellSelector)>0)/sum(cellSelector);
end

ensAve_noStim = nanmean(grayCtl_noStim);
ensSTD_noStim = nanstd(grayCtl_noStim)/sqrt(length(grayCtl_noStim));
ensAve_vis = nanmean(grayCtl_vis);
ensSTD_vis = nanstd(grayCtl_vis)/sqrt(length(grayCtl_vis));

figure(324); clf;
subplot(1,2,1); hold on

plot(0.95+0.1*rand(length(grayCtl_vis),1),grayCtl_vis,'.','markersize',16)
errorbar(1,ensAve_vis,ensSTD_vis,'linewidth',1.5,'color','k')

plot(1.25+0.1*rand(length(grayCtl_noStim),1),grayCtl_noStim,'.','markersize',16)
errorbar(1.3,ensAve_noStim,ensSTD_noStim,'linewidth',1.5,'color','k')
plot([0.5 1.5],0*[0.5 1.5],'k--','markersize',16)
xlim([0.8 1.5])
set(gca,'fontsize',16)
xticks([1 1.3])
xticklabels({'Gray (vis epochs)','Gray (no stim)'})

subplot(1,2,2); hold on
for ii = 1:18
    plot([1 2],[grayCtl_vis(ii) grayCtl_noStim(ii)],'k.-','linewidth',1.5,'markersize',16)
end
boxplot([grayCtl_vis grayCtl_noStim],[ones(length(grayCtl_vis),1); 2*ones(length(grayCtl_noStim),1)])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Gray screen (vis epochs)','Gray (no stim)'})
set(gca,'fontsize',16)

% Fig3_ctl(cellTable,cellCond)


%%
figure(); clf; hold on;

for ii = 1:18
    plot([1 2],[perPost_noStim(ii) perPost_vis(ii)],'k.-','linewidth',1.5,'markersize',16)
end
boxplot([perPost_noStim perPost_vis],[ones(length(perPost_noStim),1); 2*ones(length(perPost_vis),1)])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Gray (vis)','Gray (no stim)'})
set(gca,'fontsize',16)
ylabel('Fraction Act.')

