%%
% Plots the percent of cells activated/suppressed
%
% Run cellByCellAnalysis_GH to use this function
%%
function FigPercentAct(cellTable,cellCond,thresVal)

distBins = [15:15:250];
plotDist = distBins(1:end-1) + 15/2;

cellAct=zeros(length(distBins)-1,1);
cellSupp=zeros(length(distBins)-1,1);
for ll = 1:length(distBins)-1
    cellSelectorDist = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
    
    cellSelector = cellSelectorDist & cellCond;
    
    cellAct(ll) = sum(cellTable.dff(cellSelector)>thresVal)/sum(cellSelector);
    cellSupp(ll) = sum(cellTable.dff(cellSelector)<-thresVal)/sum(cellSelector);
end

figure(); clf; hold on;
plot(plotDist,cellAct*100,'linewidth',1.5)
plot(plotDist,cellSupp*100,'linewidth',1.5)
if thresVal == 0
    plot(plotDist,50+0*plotDist,'k--','linewidth',1.5)
end
set(gca,'fontsize',16)
legend('Activated','Suppressed')
xlabel('Min Dist to Stim Location')

end

