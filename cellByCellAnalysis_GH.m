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
           (All.out.exp.visID == 1 |  All.out.exp.visID == 0 ); % not sure what this is   
        %ismember(All.out.exp.visID, [0 1]);
                
        tempTrial = tempTrial + sum(trialsToUse);
        %% Cells to use
        holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
        
        % tg = All.out.exp.rois{holo}; % these ROIs for cells in ensemble
        % These are the cell indices for those in ensemble
        tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));
        
        offTargetRisks = All.out.anal.offTargetRisk(holo,:);
        ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
        pVis = All.out.anal.pVisR;
        % visCells = All.out.anal.pVisR < 0.05; % visually responsive
        % cellsToUse = ~offTargetRisks & ~ROIinArtifact' & visCells;
        cellsToUse = ~ROIinArtifact' & ~All.out.anal.cellsToExclude;
                    
        tgCell = zeros(length(cellsToUse),1);
        tgCell(tg_adj) = true;
        numOfTgs(ensNum) = length(tg_adj);
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
%         if norm(dffCellResp-dffCellResp2)~=0
%            plot(dffCellResp,dffCellResp2)
%            'here';
%         end
        
        
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
        
        tempIndices = find(dffCellResp>2 & cellDist<20 & ~offTargetRisks(cellsToUse)');
        if ~isempty(tempIndices) && 0 
            tempIndices %(664)
            %%
%             tempIndices = tg_adj();
            tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
            
            recWinSec=opts.recWinRange+All.out.exp.visStart;
            winToUse = round(recWinSec*All.out.info.FR);
            bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
            baselineEst = squeeze(mean(tracesHolodfData(tempIndices(1),bwinToUse(1):bwinToUse(2),:)));
            
            figure(131432499); clf; 
            subplot(2,length(tg_adj),1); hold on;
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
            ylim([-1 5])
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
            tgCell(cellsToUse) ensNum*ones(sum(cellsToUse),1) expNum*ones(sum(cellsToUse),1) cellEnsPO];
        ensNum = ensNum + 1;
    end
    
    if ~isempty(cellsToUse)
        expNum = expNum +1;
    end
    
    fprintf('Exp %d out of %d done loading\n',outer_loop, length(loadList_all));
end

%% Turn matrix into a table (easier to access specific columns by name)

cellTable = array2table(cellData,'VariableNames',{'dff','cellDist','cellPO','cellOSI','cellEnsMaxD',...
    'cellMeanEnsOSI','cellEnsOSI','visP','offTarget','tgCell','ensNum','expNum','ensPO'});

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

%% Figure Space: Effect of the spread of ensemble
ensResp = FigSpace(cellTable);

%% Figure 5: Iso vs. ortho

Fig5(cellTable,cellCondTuned);
Fig5_cbc(cellTable,cellCondTuned);

%% Figure 6: Tight co-tuned investigation
Fig6(cellTable,cellCond)
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
