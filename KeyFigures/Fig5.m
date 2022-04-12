function [cellOrisDiff] = Fig5(cellTable)

%%

oriVals = [NaN 0:45:315];
cellOris = oriVals(cellTable.cellPO)';

cellOrisDiff = abs(cellOris-cellTable.ensPO);
cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
cellOrisDiff(cellOrisDiff==135)=45;
cellOrisDiff(cellOrisDiff==180)=0;

%%
distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

cellDistDataAveOrtho=[];
cellDistDataAveIso=[];
numEnsPref = zeros(length(distBins)-1,1);
numEnsOrtho = zeros(length(distBins)-1,1);
for ii = 1:160
    for ll = 1:length(distBins)-1
        cellSelectorConds = cellTable.ensNum == ii & ...
            cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1) ...
            & cellTable.offTarget==0;
        cellSelectorEns = cellTable.cellEnsOSI>0.7 & cellTable.cellMeanEnsOSI>0.5;
        cellSelectorOriIso = cellOrisDiff == 0 & cellTable.cellOSI>0.25 & cellTable.visP<0.05;
        cellSelectorOriOrtho = cellOrisDiff == 90 & cellTable.cellOSI>0.25 & cellTable.visP<0.05;
        
        cellSelectorIso= cellSelectorConds & cellSelectorEns & cellSelectorOriIso;
        cellSelectorOrtho= cellSelectorConds & cellSelectorEns & cellSelectorOriOrtho;
       
        cellDistDataAveIso(ll,ii) = nanmean(cellTable.dff(cellSelectorIso));
        cellDistDataAveOrtho(ll,ii) = nanmean(cellTable.dff(cellSelectorOrtho));
        
        if sum(cellSelectorIso)~=0
            numEnsPref(ll) = numEnsPref(ll)+1;
        end
        
        if sum(cellSelectorOrtho)~=0
            numEnsOrtho(ll) = numEnsOrtho(ll)+1;
        end
    end
end

respAveIso = nanmean(cellDistDataAveIso,2);
respStdErrIso = nanstd(cellDistDataAveIso,[],2)./sqrt(numEnsPref);
respAveOrtho = nanmean(cellDistDataAveOrtho,2);
respStdErrOrtho = nanstd(cellDistDataAveOrtho,[],2)./sqrt(numEnsOrtho);

figure(); clf; hold on;
hLeg(1) = plot(plotDist,respAveIso,'.-','linewidth',1.5,'markersize',15,'color',[0 0.447 0.741]);
errorbar(plotDist,respAveIso,respStdErrIso,'linewidth',1.5,'color',[0 0.447 0.741])
hLeg(2) = plot(plotDist,respAveOrtho,'.-','linewidth',1.5,'markersize',15,'color',[0.494 0.184 0.556]);
errorbar(plotDist,respAveOrtho,respStdErrOrtho,'linewidth',1.5,'color',[0.494 0.184 0.556])
plot([0 250],0*[0 250],'k--')
set(gca,'fontsize',16)
xlim([0 150])
legend(hLeg,{'Iso','Ortho'})

end


