%%
% Plots the percent of cells activated/suppressed
%
% Run cellByCellAnalysis_GH to use this function
%%
function SFigPercentActEns(cellTable,cellCond,thresVal)


totalNumEns = cellTable.ensNum(end);
distBins = [15:15:250];
plotDist = distBins(1:end-1) + 15/2;

cellAct=zeros(length(distBins)-1,1);
cellSupp=zeros(length(distBins)-1,1);

% Loop over all ensembles
for ii = 1:totalNumEns
    for ll = 1:length(distBins)-1
        cellSelectorDist = cellTable.ensNum == ii &...
            cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
        
        cellSelector = cellSelectorDist & cellCond;
        
        cellAct(ll,ii) = sum(cellTable.dff(cellSelector)>thresVal)/sum(~isnan(cellTable.dff(cellSelector)));
        cellSupp(ll,ii) = sum(cellTable.dff(cellSelector)<-thresVal)/sum(~isnan(cellTable.dff(cellSelector)));
    end
end

%%
% Average across ensembles
actAve = nanmean(cellAct,2);
suppAve = nanmean(cellSupp,2);
actStdErr = nanstd(cellAct,[],2)/sqrt(size(cellAct,2));
suppStdErr = nanstd(cellSupp,[],2)/sqrt(size(cellSupp,2));

%%

figure(); clf; 
subplot(1,2,1); hold on;
% plot(plotDist,actAve*100,'linewidth',1.5,'color',[0.8500    0.3250    0.0980])
errorbar(plotDist,actAve*100,actStdErr*100,'color',[0.8500    0.3250    0.0980],'linewidth',2,'CapSize',0);
% plot(plotDist,suppAve*100,'linewidth',1.5,'color',[0    0.4470    0.7410])
errorbar(plotDist,suppAve*100,suppStdErr*100,'color',[0    0.4470    0.7410],'linewidth',2,'CapSize',0);
if thresVal == 0
    plot([0 250],50+0*[0 250],'k--','linewidth',1.5)
end
xlim([0 250])
set(gca,'fontsize',16)
legend('Activated','Suppressed')
xlabel('Min Dist to Stim Location')
ylabel('Percent of cells')

subplot(1,2,2)
histogram(cellAct(1,:)*100)
set(gca,'fontsize',16)
xlabel('Percent of cells activated in first bin')
xticks([0:25:100])
ylabel('Number of Ensembles')

end

