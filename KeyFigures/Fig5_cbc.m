function [cellResponsesIso,cellResponsesOrtho] = Fig5_cbc(cellTable,cellOrisDiff)

distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

cellDistDataAveOrtho=[];
cellDistDataAvePref=[];
numEnsIso = zeros(length(distBins)-1,1);
numEnsOrtho = zeros(length(distBins)-1,1);

respAveIso=[]; respStdErrIso =[]; respMedIso = [];
respAveOrtho=[]; respStdErrOrtho = []; respMedOrtho = [];

for ll = 1:length(distBins)-1
    cellSelectorConds = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1) ...
        & cellTable.offTarget==0;
    cellSelectorEns = cellTable.cellEnsOSI>0.7 & cellTable.cellMeanEnsOSI>0.5;
    cellSelectorOriIso = cellOrisDiff == 0 & cellTable.cellOSI>0.25 & cellTable.visP<0.05;
    cellSelectorOriOrtho = cellOrisDiff == 90 & cellTable.cellOSI>0.25 & cellTable.visP<0.05;
    
    cellSelectorIso= cellSelectorConds & cellSelectorEns & cellSelectorOriIso;
    cellSelectorOrtho= cellSelectorConds & cellSelectorEns & cellSelectorOriOrtho;
    
    respAveIso(ll) = nanmean(cellTable.dff(cellSelectorIso));
    respStdErrIso(ll) = nanstd(cellTable.dff(cellSelectorIso))./sqrt(sum(cellSelectorIso));
    respMedIso(ll) = nanmedian(cellTable.dff(cellSelectorIso));
    
    respAveOrtho(ll) = nanmean(cellTable.dff(cellSelectorOrtho));
    respStdErrOrtho(ll) = nanstd(cellTable.dff(cellSelectorOrtho))./sqrt(sum(cellSelectorOrtho));
    respMedOrtho(ll) = nanmedian(cellTable.dff(cellSelectorOrtho));
    
    if ll == 1
        cellResponsesIso = cellTable.dff(cellSelectorIso);
        cellResponsesOrtho = cellTable.dff(cellSelectorOrtho);
    end
end

figure(23123); hold on;
hLeg(1) = plot(plotDist,respAveIso,'.-','linewidth',1.5,'markersize',15,'color',[0 0.447 0.741]);
errorbar(plotDist,respAveIso,respStdErrIso,'linewidth',1.5,'color',[0 0.447 0.741])
plot(plotDist,respMedIso,'*','linewidth',1.5,'markersize',15,'color',[0 0.447 0.741]);
hLeg(2) = plot(plotDist,respAveOrtho,'.-','linewidth',1.5,'markersize',15,'color',[0.494 0.184 0.556]);
errorbar(plotDist,respAveOrtho,respStdErrOrtho,'linewidth',1.5,'color',[0.494 0.184 0.556])
plot(plotDist,respMedOrtho,'*','linewidth',1.5,'markersize',15,'color',[0.494 0.184 0.556]);
plot([0 250],0*[0 250],'k--')
set(gca,'fontsize',16)
xlim([0 150])
legend(hLeg,{'Iso','Ortho'})


end

