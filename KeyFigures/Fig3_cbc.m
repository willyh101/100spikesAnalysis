function [outputArg1,outputArg2] = Fig3_cbc(cellTable)
% distBins = [15:15:250];
distBins = [15:15:250];
plotDist = distBins(1:end-1) + 15/2;

cellDistDataAve = zeros(1,length(distBins)-1);
cellDistDataMedian = zeros(1,length(distBins)-1);
cellDistDataErr = zeros(1,length(distBins)-1);

cellDistData=[];
for ll = 1:length(distBins)-1
    cellSelector = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1) ... 
        & cellTable.offTarget==0;
    
    cellDistDataAve(ll) = nanmean(cellTable.dff(cellSelector));
    cellDistDataMedian(ll) = nanmedian(cellTable.dff(cellSelector));
    cellDistDataErr(ll) = nanstd(cellTable.dff(cellSelector))/sqrt(sum(cellSelector));
end


figure(); clf;
subplot(1,2,1);
plot(cellTable.cellDist(cellTable.offTarget==0),cellTable.dff(cellTable.offTarget==0),'.')
set(gca,'fontsize',16)
xlabel('Min Dist to Stim Location')
xlim([15 50])
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
legend(leg,{'Mean','Median'})

end

