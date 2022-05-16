function Fig3_noStimCtl(cellTable)

%%
binSize = 15;
distBins = [0:binSize:200];
plotDist = distBins(1:end-1) + binSize/2;
cellDistDataAve=zeros(length(distBins)-1,cellTable.expNum(end));
for ii = 1:cellTable.expNum(end)
    
    expSelector = cellTable.expNum ==ii;
    posEns = unique(cellTable.ensNum(expSelector));
    
    cellDistDataAve1 = zeros(length(distBins)-1,length(posEns));
    for ee = 1:length(posEns)
        ensSelector = cellTable.ensNum == posEns(ee);
        for ll = 1:length(distBins)-1
            cellSelectorDist = cellTable.cellDist>=distBins(ll) & cellTable.cellDist<distBins(ll+1);
            cellSelector = cellSelectorDist & ensSelector;
             
            cellDistDataAve1(ll,ee) = nanmean(cellTable.('ctlResp_noStim')(cellSelector));
        end
    end
    cellDistDataAve(:,ii) = nanmean(cellDistDataAve1,2);
end
respAve = nanmean(cellDistDataAve,2);
respStdErr = nanstd(cellDistDataAve,[],2)/sqrt(size(cellDistDataAve,2));

figure(4253); clf; hold on;
hLeg(1) = plot(plotDist,respAve,'.-','color',[0 0.4470 0.7410],'linewidth',1.5,'markersize',15);
errorbar(plotDist,respAve,respStdErr,'color',[0 0.4470 0.7410],'linewidth',1.5)

plot([0 200],0*[0 200],'k--')
set(gca,'fontsize',16)
xlabel('Min dist to (fake) stim')
ylabel('Mean \DeltaF/F')
ylim([-0.05 0.02])
title('All cells')

legend(hLeg,'Gray screen (no stim)')
end

