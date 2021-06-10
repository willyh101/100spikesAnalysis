colorListOri = colorMapPicker(numel(diffsPossible),'hsv');

for k = 1:numel(diffsPossible)
    f = figure(110);
    clf
    
    dat = popToPlot(:,:,k);
    meanDat = nanmean(dat);
    stdDat = nanstd(dat);
    numpDat = sum(~isnan(dat));
    semDat = stdDat./sqrt(numpDat);
    distBinSize = distBins(2)-distBins(1);
    e = errorbar(distBins(2:end)-distBinSize/2,meanDat,semDat,'linewidth',2,'color',colorListOri{k});
    xlim([0 opts.distBins(end)])
    ylim([-0.1 0.3])
    title(['Cells Preferred Angle \Delta' num2str(diffsPossible(k)) '\circ'])
    ylabel('Mean z-scored dF/F')
    xlabel('Distance from Ensemble (\mum)')
    yline(0,'--k','LineWidth',1)
    box off
    
    fname = ['cells pref angle delta ori ' num2str(k)];
    saveas(f, ['c:/users/will/desktop/cortex club/' fname '.png'])
    
%     pause
end



