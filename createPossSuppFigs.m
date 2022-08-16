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

%%
load('./compressedData/cellTable220816.mat')

%% Cell conditions used in the functions
cellCond = cellTable.offTarget==0; 
cellCondTuned = cellTable.offTarget==0 & cellTable.visP<0.05 & cellTable.cellOSI > 0.25;
cellCondNonVis = cellTable.offTarget==0 & cellTable.visP>0.05;

%% Possible supplemental figure: min dist vs. dF/F plot for gray screen
SFig_noStimCtl(cellTable)

%% Supplemental figure: min dist vs. dF/F plot for different bin widths
fprintf('Creating min dist vs. dF/F plot fits for different bins\n')
SFig_dataFit(cellTable,cellCond)
% SFig_cbc_dataFit(cellTable,cellCond)

%% Percent of cells suppressed vs. activated across all ensembles

totalNumEns = cellTable.ensNum(end);
allResp =zeros(totalNumEns,1);
for ii = 1:totalNumEns
    cellSelector = cellCond & cellTable.ensNum == ii & cellTable.cellDist<inf;
    allResp(ii) = sum(cellTable.dff(cellSelector)<0)/sum(~isnan(cellTable.dff(cellSelector)));
end

figure(7773); clf; 
subplot(1,2,1); hold on
violins = violinplot(allResp*100);
plot([0 2],50+[0 0],'k--','linewidth',1.5)
ylabel('% of cells suppressed')
set(gca,'fontsize',16)

subplot(1,2,2); hold on
histogram(allResp*100)
xlabel('% of cells suppressed')
ylabel('Number of Ensembles')
set(gca,'fontsize',16)
%%
cellSelector = cellCond & cellTable.cellDist<200;
allResp = cellTable.dff(cellSelector);
allResp(isnan(allResp))=[];
temp1=sum(allResp<0)/sum(~isnan(allResp));
temp2=sum(allResp>0)/sum(~isnan(allResp));

fprintf('Percent of cells activated: %.2f\n',temp2*100)
figure(7774); clf; hold on
bar([temp1;temp2]*100)
set(gca,'fontsize',16)
xticks([1 2])
xticklabels({'Suppressed','Activated'})
ylabel('Percent of cells')
plot([0 3],50+[0 0],'k--','linewidth',1.5)
xlim([0 3])

%%

figure();

subplot(1,2,1)
plot(unique(cellTable.cellEnsMeaD,'stable'),unique(cellTable.cellEnsOSI,'stable'),'.','markersize',16)
set(gca,'fontsize',16)
xlabel('Ens Spread')
ylabel('Ens OSI')
ylim([0 1])

subplot(1,2,2)
plot(unique(cellTable.cellEnsMeaD,'stable'),unique(cellTable.cellMeanEnsOSI,'stable'),'.','markersize',16)
set(gca,'fontsize',16)
xlabel('Ens Spread')
ylabel('Mean Ens OSI')
ylim([0 1])

%%

%% Figure plotting the percent activated/suppressed as a function of dist
% Third input is threshold value
SFigPercentAct(cellTable,cellCond,0)

%%
SFigPercentActEns(cellTable,cellCond,0)

%% Target cell statistics
SFigTargets(cellTable)

%% Figure illustrating the total number of cells suppressed
% clc;
fprintf('-------------------------\n')
SFigCellByCell(cellTable,cellCond)

%%

totalExps = cellTable.expNum(end);
grayCtl_noStim = [];
perPost_noStim = [];
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

