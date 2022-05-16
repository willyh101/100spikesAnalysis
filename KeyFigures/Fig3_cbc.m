%%
% Reproduces the minimal distance to the ensemble plot by not averaging
% across all ensembles first
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig3_cbc(cellTable,cellCond)

distBins = [15:15:250];
plotDist = distBins(1:end-1) + 15/2;

cellDistDataAve = zeros(1,length(distBins)-1);
cellDistDataMedian = zeros(1,length(distBins)-1);
cellDistDataErr = zeros(1,length(distBins)-1);
% Loop over all distances
for ll = 1:length(distBins)-1
    cellSelectorDist = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
    cellSelector = cellSelectorDist & cellCond;
    
    % Average across distances
    cellDistDataAve(ll) = nanmean(cellTable.dff(cellSelector));
    cellDistDataMedian(ll) = nanmedian(cellTable.dff(cellSelector));
    cellDistDataErr(ll) = nanstd(cellTable.dff(cellSelector))/sqrt(sum(cellSelector));
end

% Plot the results
figure(); clf;
subplot(1,2,1);
plot(cellTable.cellDist(cellCond),cellTable.dff(cellCond),'.')
set(gca,'fontsize',16)
xlabel('Min Dist to Stim Location')
xlim([0 250])
ylabel('\DeltaF/F')

subplot(1,2,2); hold on;
leg(1) = plot(plotDist,cellDistDataAve,'k.-','markersize',16);
errorbar(plotDist,cellDistDataAve,cellDistDataErr,'k','linewidth',1.5)
leg(2) = plot(plotDist,cellDistDataMedian,'*','markersize',16);
plot([0 distBins(end)],0*[0 distBins(end)],'k--');
set(gca,'fontsize',16)
xlabel('Min Dist to Stim Location')
ylabel('Mean Evoked \DeltaF/F')
xlim([0 250])
xticks([0:25:250])
legend(leg,{'Mean','Median'})

end

