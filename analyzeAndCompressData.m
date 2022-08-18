%%
% Analyzes the experimental data to produce cellTable, which has the core
% statistics we are interested in (e.g., dF/F, OSI)
%%
clear; close all; clc;

%%
% Adds all subfolders to the path
restoredefaultpath;
folder = fileparts(which('cellByCellAnalysis_GH.m')); 
addpath(genpath(folder));
rmpath(folder);

%% Specify data location and loadList

% loadPath = '/Users/gregoryhandy/Research_Local/outfiles/'; % where ever your files are
loadPath = '/Users/greghandy/Research_Local_v2/'; % where ever your files are

% Options: short (for debugging), all (all avaliable outfiles), used (data
% files currently being used here; created for speed)
loadList_all = oriLoadList_GH('used');

%% Loop over each experiment
ensNum = 1; expNum = 1;
cellData = []; expData=[];
cellID=0;
%%
for outer_loop = 1:length(loadList_all) %10 for 60 offTarget
    
    %% Load the data from a specific experiment (code from Will)
    loadList = loadList_all(outer_loop);
    [All,opts,outVars] = loadData_GH(loadList,loadPath);
    
    %% Mouse name
    mouseName = loadList_all(outer_loop);
    mouseName = mouseName{1};
    tempIndices = find(mouseName=='_');
    mouseName([1:tempIndices(1) tempIndices(end):end])='';

    %% Get orientation tuning and OSI from dataset
    [All, outVars] = getTuningCurve(All, opts, outVars);
    [All, outVars] = calcOSI(All, outVars);
    [All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
    [All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically
    
    [All, outVars] = defineDistanceTypes(All, outVars);
    All.out.exp.ensMaxD =  outVars.ensMaxD;
    All.out.exp.ensMeaD =  outVars.ensMeaD;
           
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
        cellsToUse = ~All.out.anal.cellsToExclude & ~ROIinArtifact';% & ;
        
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
                
        % Method 2: use results from anal 
        tempResp = squeeze(All.out.anal.respMat(s,cellsToUse));
        tempBase = squeeze(All.out.anal.baseMat(s,cellsToUse));
        dffCellResp2 = (tempResp-tempBase)';
        
        % Methods should be identical
        if norm(dffCellResp(~isnan(dffCellResp))-dffCellResp2(~isnan(dffCellResp)))~=0 ...
                || sum(isnan(dffCellResp))~=sum(isnan(dffCellResp2))
            plot(dffCellResp,dffCellResp2)
            error('Something went wrong!');
        end
                
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
                
        % The gray screen should be the first stimID on this list
        % Here, we do a quick check for this
        if unique(All.out.anal.minDistbyHolo(1,:)) ~=0 || All.out.exp.stimParams.roi{1}~=0
            error('Something wrong with minDistbyHolo');
        end
        
        %% Careful indexing into All.out.anal.minDistbyHolo
        cellDist = All.out.anal.minDistbyHolo(ensID,cellsToUse)';
        
        %% These are not ensemble-dependent
        cellPO = All.out.anal.prefOri(cellsToUse)';
        cellOSI = All.out.anal.osi(cellsToUse)';
        
        %% Shared across all cells 
        cellEnsMaxD = All.out.exp.ensMaxD(holo)*ones(sum(cellsToUse),1);
        cellEnsMeaD = All.out.exp.ensMeaD(holo)*ones(sum(cellsToUse),1);
        cellMeanEnsOSI = All.out.exp.meanEnsOSI(holo)*ones(sum(cellsToUse),1);
        cellEnsOSI = All.out.exp.ensOSI(holo)*ones(sum(cellsToUse),1);
        cellEnsPO = All.out.exp.ensPO(holo)*ones(sum(cellsToUse),1);
                
        %% Store the data
        if ii==1
            IDStart = cellID(end)+1;
            cellID = (IDStart:(IDStart-1)+sum(cellsToUse))';
        end
               
        cellData = [cellData; dffCellResp cellDist cellPO cellOSI ...
            cellEnsMaxD cellEnsMeaD cellMeanEnsOSI cellEnsOSI pVis(cellsToUse)' offTargetRisks(cellsToUse)'...
            tgCell(cellsToUse) ensNum*ones(sum(cellsToUse),1) expNum*ones(sum(cellsToUse),1) cellEnsPO,...
            dffCellCtlResp_noStim cellID];
        mouseNames{ensNum} = mouseName;
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
    'cellEnsMeaD', 'cellMeanEnsOSI','cellEnsOSI','visP','offTarget','tgCell','ensNum','expNum','ensPO',...
    'ctlResp_noStim','cellID'});

%% Iso vs. ortho calculation
oriVals = [NaN 0:45:315];
cellOris = oriVals(cellTable.cellPO)';
cellOrisDiff = abs(cellOris-cellTable.ensPO);
cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
cellOrisDiff(cellOrisDiff==135)=45;
cellOrisDiff(cellOrisDiff==180)=0;
cellTable.cellOrisDiff = cellOrisDiff;

%% Optional: After the code runs, save cellTable to the compressedData folder 

% clearvars -except cellTable mouseNames
% save('cellTableTodaysDate')

