%%
% Reproduces the minimal distance to the ensemble plot
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%%
function Fig3(cellTable,cellCond)

totalNumEns = cellTable.ensNum(end);

distBins = [15:15:250];
plotDist = distBins(1:end-1) + 15/2;

cellDistDataAve=zeros(length(distBins)-1,totalNumEns);
% Loop over all ensembles
for ii = 1:totalNumEns
    % Loop over all distances
    for ll = 1:length(distBins)-1
        cellSelectorDist = cellTable.ensNum == ii & ...
            cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
        
        cellSelector = cellSelectorDist & cellCond;
        
        cellDistDataAve(ll,ii) = nanmean(cellTable.('dff')(cellSelector));
    end
end
% Average across ensembles
respAve = nanmean(cellDistDataAve,2);
respStdErr = nanstd(cellDistDataAve,[],2)/sqrt(size(cellDistDataAve,2));

% Plot the results
figure(); clf; hold on;
plot(plotDist,respAve,'k.-','linewidth',1.5,'markersize',15)
errorbar(plotDist,respAve,respStdErr,'k','linewidth',1.5)
plot([0 250],0*[0 250],'k--')
xlim([0 250])
maxVal = abs(max(respAve+respStdErr));
minVal = abs(min(respAve-respStdErr));
ylim([-round(max(minVal,maxVal),2) round(max(minVal,maxVal),2)])
set(gca,'fontsize',16)
ylabel('Mean evoked \DeltaF/F')
xlabel('Min Dist to Stim Loc')

end

