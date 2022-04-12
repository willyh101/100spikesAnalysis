function [outputArg1,outputArg2] = Fig3(cellTable)

distBins = [15:15:250];
plotDist = distBins(1:end-1) + 15/2;

cellDistDataAve=[];
for ii = 1:160
    for ll = 1:length(distBins)-1
        cellSelector = cellTable.ensNum == ii & ...
            cellTable.('cellDist')>distBins(ll) & cellTable.('cellDist')<distBins(ll+1) ...
            & cellTable.('offTarget')==0;
        
        cellDistDataAve(ll,ii) = nanmean(cellTable.('dff')(cellSelector));
    end
end
respAve = nanmean(cellDistDataAve,2);
respStdErr = nanstd(cellDistDataAve,[],2)/sqrt(size(cellDistDataAve,2));

figure(); clf; hold on;
plot(plotDist,respAve,'k.-','linewidth',1.5,'markersize',15)
errorbar(plotDist,respAve,respStdErr,'k','linewidth',1.5)
plot([0 250],0*[0 250],'k--')
set(gca,'fontsize',16)
xlim([0 250])
ylim([-0.05 0.05])

end

