%%
%
%
% Run cellByCellAnalysis_GH to use this function
%%
function FigCellByCell(cellTable, cellCond)

thresVal = [0 0.05];

for thresLoop = 1:2
    suppFrac = zeros(cellTable.ensNum(end),1);
    actFrac = zeros(cellTable.ensNum(end),1);
    for ii = 1:cellTable.ensNum(end)
        cellSelectorAll = cellCond & cellTable.ensNum == ii & ~isnan(cellTable.dff);
        cellRespTemp = cellTable.dff(cellSelectorAll);
        
        actFrac(ii) = sum(cellRespTemp>thresVal(thresLoop))/sum(cellSelectorAll);
        suppFrac(ii) = sum(cellRespTemp<-thresVal(thresLoop))/sum(cellSelectorAll);
    end
    
    if thresLoop == 1
        figure();clf;
        subplot(1,2,1); hold on;
        histogram(actFrac*100,[20:5:70])
        plot(mean(actFrac*100),[60],'*','markersize',15)
        plot(50+[0 0],[0 60],'k--','linewidth',1.5)
        set(gca,'fontsize',16)
        xlabel('% of Population Activated')
        title(sprintf('Threshold: %.2f',thresVal(thresLoop)))
        
        subplot(1,2,2); hold on;
        boxplot(suppFrac*100)
        set(gca,'fontsize',16)
        ylabel('% of Population Suppressed')
        xticklabels('')
        
    else
        figure(); clf;
        subplot(1,2,1); hold on;
        histogram((suppFrac-actFrac)*100)
        plot(mean((suppFrac-actFrac)*100),50,'*','markersize',15)
        plot(0+[0 0],[0 50],'k--','linewidth',1.5)
        set(gca,'fontsize',16)
        xlabel('Suppressed % - Activated %')
        title(sprintf('Threshold: %.2f',thresVal(thresLoop)))
        
        subplot(1,2,2); hold on;
        boxplot((suppFrac-actFrac)*100)
        set(gca,'fontsize',16)
        ylabel('Suppressed % - Activated %')
        xticklabels('')
    end
    
    fprintf('Thres: %0.2f \n',thresVal)
    fprintf('Percent of cells suppressed: %.2f%% %s %.2f%%\n',mean(suppFrac)*100,char(177),std(suppFrac)/sqrt(cellTable.ensNum(end))*100)
    fprintf('Percent of cells activated: %.2f%% %s %.2f%%\n',mean(actFrac)*100,char(177),std(actFrac)/sqrt(cellTable.ensNum(end))*100)
    
end

end

