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

% loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/'; % where ever your files are
loadPath = '/Users/greghandy/Research_Local_v2/'; % where ever your files are
% loadList_all = oriLoadList_GH('all');
loadList_all = oriLoadList_GH('used');

ensNum = 1;
expNum = 1;
cellData = [];
%%
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
            All.out.exp.stimSuccessTrial &... %All.out.exp.visID == 1;
           (All.out.exp.visID == 1 |  All.out.exp.visID == 0 ); % not sure what this is   
        %ismember(All.out.exp.visID, [0 1]);
                
        tempTrial = tempTrial + sum(trialsToUse);
        %% Cells to use
        holo = All.out.exp.stimParams.roi{s};  % tracks stimID -> hologram
        
%         tg = All.out.exp.rois{holo}; % these are cells in ensemble
        % fixed finding the cells in the ensembles? 
        % Not used for this analysis
        tg_adj = All.out.exp.holoTargets{holo}(~isnan(All.out.exp.holoTargets{holo}));
        
        offTargetRisks = All.out.anal.offTargetRisk(holo,:);
        ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
        pVis = All.out.anal.pVisR;
%         visCells = All.out.anal.pVisR < 0.05; % visually responsive
%         cellsToUse = ~offTargetRisks & ~ROIinArtifact' & visCells;
        cellsToUse = ~ROIinArtifact' & ~All.out.anal.cellsToExclude;
                    
        
        tgCell = zeros(length(cellsToUse),1);
        tgCell(tg_adj) = 1;
        %% Baseline the df data
        % Options: dfData, zdfData, allData
        
        % Method 1: 
        %recWinSec=opts.recWinRange+All.out.exp.visStart;
        %winToUse = round(recWinSec*All.out.info.FR);
        %bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
        %tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
        %dfCellResp = mean((mean(tracesHolodfData(:,winToUse(1):winToUse(end),:),2)),3);
        %baselineEst = mean((mean(tracesHolodfData(:,bwinToUse(1):bwinToUse(2),:),2)),3);
        %dffCellResp = (dfCellResp-baselineEst);
         
        % Method 2: use results from anal 
        tempResp = squeeze(All.out.anal.respMat(s,cellsToUse));
        tempBase = squeeze(All.out.anal.baseMat(s,cellsToUse));
        dffCellResp = (tempResp-tempBase)';
        
        % Methods should be identical, but something is off...
        %if norm((temp-temp_base)-dffCellResp')~=0
           %plot(temp-temp_base,dffCellResp)
           %'here';
        %end
        
        
        %% Unclear how to index into All.out.anal.minDistbyHolo
        
        % All.out.exp.ensMaxD has 20 entires but All.out.anal.minDistbyHolo
        % has 21 entries
        % All.out.anal.minDistbyHolo(1,:) is all 0's so using holo seems
        % wrong. So instead, use ensID (or s) from above
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
%         if ~isempty(tempIndices)
%             tempIndices %(664)
%             %%
% %             tempIndices = tg_adj();
%             tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
%             
%             recWinSec=opts.recWinRange+All.out.exp.visStart;
%             winToUse = round(recWinSec*All.out.info.FR);
%             bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
%             baselineEst = squeeze(mean(tracesHolodfData(tempIndices(1),bwinToUse(1):bwinToUse(2),:)));
%             
%             figure(1); clf; hold on;
%             baselinedResp = zeros(12,24);
%             for gg = 1:size(tracesHolodfData,3)
%                 baselinedResp(gg,:) = tracesHolodfData(tempIndices(1),:,gg)-baselineEst(gg);
%                 plot(baselinedResp(gg,:),'color',[0 0 0 0.3])
%             end
%             
%             plot(mean(baselinedResp),'k','linewidth',2)
%             plot(round(All.out.exp.visStart*All.out.info.FR)+[0 0],[-1 7],'k--')
%             plot(round((All.out.exp.visStart+1)*All.out.info.FR)+[0 0],[-1 7],'k--')
%           
%             tempResp = mean(baselinedResp);
%             mean(tempResp(winToUse(1):winToUse(2)))
%             dffCellResp(tempIndices(1))
%             %%
%         end
        %%
%         tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
%         recWinSec=opts.recWinRange+All.out.exp.visStart;
%         winToUse = round(recWinSec*All.out.info.FR);
%         bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
%         figure(4); clf; hold on;
%         for stimLoop = 1:length(tg_adj)
%             
%             baselineEst = squeeze(mean(tracesHolodfData(tg_adj(stimLoop),bwinToUse(1):bwinToUse(2),:)));
%             subplot(1,7,stimLoop); hold on;
%             for gg = 1:size(tracesHolodfData,3)
%                 baselinedResp(gg,:) = tracesHolodfData(tg_adj(stimLoop),:,gg)-baselineEst(gg);
%                 plot(baselinedResp(gg,:),'color',[0 0 0 0.3])
%             end
%             plot(mean(baselinedResp),'k','linewidth',2)
%         end
        
        %%
        
        
        cellData = [cellData; dffCellResp...
            cellDist cellPO cellOSI cellEnsMaxD cellMeanEnsOSI cellEnsOSI pVis(cellsToUse)' offTargetRisks(cellsToUse)'...
            tgCell(cellsToUse) ensNum*ones(sum(cellsToUse),1) expNum*ones(sum(cellsToUse),1) cellEnsPO];
        ensNum = ensNum + 1;
    end
    
    if ~isempty(cellsToUse)
        expNum = expNum +1;
    end
    
    fprintf('Exp %d out of %d done loading\n',outer_loop, length(loadList_all));
end

%%

cellTable = array2table(cellData,'VariableNames',{'dff','cellDist','cellPO','cellOSI','cellEnsMaxD',...
    'cellMeanEnsOSI','cellEnsOSI','visP','offTarget','tgCell','ensNum','expNum','ensPO'});

%% Target cell statistics
tgCellDists = [];
tgRespAve = [];
tgRespErr=[];
for ii = 1:160
    tgSelector = cellTable.ensNum == ii & cellTable.tgCell == 1;
    
    numTg(ii) = sum(tgSelector);
    tgRespAve(ii) = nanmean(cellTable.dff(tgSelector));
    tgRespErr(ii) = nanstd(cellTable.dff(tgSelector))/sqrt(sum(tgSelector));
end
    
tgSelectorAll = cellTable.tgCell == 1;
tgRespAll = cellTable.dff(tgSelectorAll);
%%

figure(22); clf;
subplot(3,2,1)
boxplot(tgRespAll)
ylabel('\DeltaF/F')
set(gca,'fontsize',16)

subplot(3,2,2)
histogram(tgRespAll)
xlabel('\DeltaF/F')
set(gca,'fontsize',16)

subplot(3,1,2)
plot(numTg,'.-')
ylabel('Number of Targets')
xlabel('Ens Number')
set(gca,'fontsize',16)

subplot(3,1,3)
errorbar(tgRespAve,tgRespErr)
set(gca,'fontsize',16)
ylabel('Mean \DeltaF/F')
xlabel('Ens Number')

%%
cellRespAll=[];
actFrac = zeros(cellTable.ensNum(end),1);
for ii = 1:cellTable.ensNum(end)
    cellSelectorAll = cellTable.offTarget==0 & cellTable.ensNum == ii ;
    cellRespTemp = cellTable.dff(cellSelectorAll);
    actFrac(ii) = sum(cellRespTemp>0)/sum(cellSelectorAll);
end
figure(11);clf; hold on;
histogram(actFrac*100,[20:5:70])
plot(50+[0 0],[0 50],'k--','linewidth',1.5)
set(gca,'fontsize',16)
xlabel('Percent of Population Activated')

clc;
fprintf('Percent of cells activated: %.2f%% %s %.2f%%\n',mean(actFrac)*100,char(177),std(actFrac)/sqrt(cellTable.ensNum(end))*100) 

%%


ensResp = FigSpace(cellTable);


figure(3423134); clf; hold on;
plot(tgRespAve,ensResp,'.','markersize',16)
plot([0 5], 0+[0 0],'k--','linewidth',1.5) 
set(gca,'fontsize',16)

fitLM2 = fit(tgRespAve',ensResp,'poly1');
plot(tgRespAve',fitLM2.p2+fitLM2.p1*tgRespAve','linewidth',2)
xlabel('Ens Response')
ylabel('Mean evoked \Delta F/F')

%% Plot the data


Fig3_cbc(cellTable)

%%

Fig3(cellTable)

%%

distBins = [15:15:250];
plotDist = distBins(1:end-1) + 15/2;
cellAct=[];
cellSupp = [];
for ll = 1:length(distBins)-1
    cellSelector = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1) ... 
        & cellTable.offTarget==0;
    
    cellAct(ll) = sum(cellTable.dff(cellSelector)>0)/sum(cellSelector);
    cellSupp(ll) = sum(cellTable.dff(cellSelector)<-0)/sum(cellSelector);
end

figure(); clf; hold on;
plot(plotDist,cellAct*100,'linewidth',1.5)
plot(plotDist,cellSupp*100,'linewidth',1.5)
plot(plotDist,50+0*plotDist,'k--','linewidth',1.5)
set(gca,'fontsize',16)
legend('Activated','Suppressed')
xlabel('Min Dist to Stim Location')
%%

cellOrisDiff = Fig5(cellTable);

%%
[cellResponsesIso, cellResponsesOrtho]=Fig5_cbc(cellTable,cellOrisDiff);

%%

figure(23434210); clf;
subplot(1,2,1); hold on;
boxplot([cellResponsesIso; cellResponsesOrtho], ...
    [ones(length(cellResponsesIso),1); 2*ones(length(cellResponsesOrtho),1)])
set(gca,'fontsize',16)

subplot(1,2,2); hold on
bar([sum(cellResponsesIso>0)/length(cellResponsesIso) ...
    sum(cellResponsesOrtho>0)/length(cellResponsesOrtho)])
plot([0 3], 0.5 + [0 0],'k--')
set(gca,'fontsize',16)




%% Tight co-tuned investigation

Fig6(cellTable,1)
Fig6IsoOrtho(cellTable,cellOrisDiff)

%%
[cellResponsesTight, cellResponsesLoose] = Fig6_cbc(cellTable,1);
%%

figure(12312423); clf;
subplot(1,2,1);
boxplot([cellResponsesTight{1}; cellResponsesLoose{1}; cellResponsesTight{2}; cellResponsesLoose{2}],...
    [ones(length(cellResponsesTight{1}),1); 2*ones(length(cellResponsesLoose{1}),1); ...
    3*ones(length(cellResponsesTight{2}),1); 4*ones(length(cellResponsesLoose{2}),1)])
set(gca,'fontsize',16)

subplot(1,2,2); hold on;
bar([sum(cellResponsesTight{1}>0)/length(cellResponsesTight{1}),...
    sum(cellResponsesLoose{1}>0)/length(cellResponsesLoose{1}),...
    sum(cellResponsesTight{2}>0)/length(cellResponsesTight{2}),...
    sum(cellResponsesLoose{2}>0)/length(cellResponsesLoose{2})])
plot([0 5],0.5+0*[0 5],'k--','linewidth',1.5)
set(gca,'fontsize',16)

