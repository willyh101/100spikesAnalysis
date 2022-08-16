%%
% Reproduces iso vs. ortho plot by not averaging across all ensembles first
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Should only included tuned and visually responsive cells
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig4_cbc(cellTable,cellCond)

distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

% Conditions for this analyis
cellSelectorEns = cellTable.cellEnsOSI>0.7 & cellTable.cellMeanEnsOSI>0.5;
cellSelectorOriIso = cellTable.cellOrisDiff == 0;
cellSelectorOriOrtho = cellTable.cellOrisDiff == 90;


respAveIso=zeros(length(distBins)-1,1); respAveOrtho=zeros(length(distBins)-1,1);
respStdErrIso=zeros(length(distBins)-1,1); respStdErrOrtho=zeros(length(distBins)-1,1);
respMedIso=zeros(length(distBins)-1,1); respMedOrtho=zeros(length(distBins)-1,1);
% Loop through all distances
for ll = 1:length(distBins)-1
    cellSelectorDist = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
    cellSelectorIso= cellSelectorDist & cellSelectorEns & cellSelectorOriIso & cellCond;
    cellSelectorOrtho= cellSelectorDist & cellSelectorEns & cellSelectorOriOrtho & cellCond;
    
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

figure(421354); hold on;
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

%% First bin analyis
figure(4213541);
subplot(1,2,1); hold on;
boxplot([cellResponsesIso; cellResponsesOrtho], ...
    [ones(length(cellResponsesIso),1); 2*ones(length(cellResponsesOrtho),1)])
set(gca,'fontsize',16)
title('First bin distributions')
ylabel('\DeltaF/F')
xticklabels({'Iso','Ortho'})
xticks([1 2])

subplot(1,2,2); hold on
bar([sum(cellResponsesIso<0)/sum(~isnan(cellResponsesIso))*100 ...
    sum(cellResponsesOrtho<0)/sum(~isnan(cellResponsesOrtho))*100])
plot([0 3], 50 + [0 0],'k--')
set(gca,'fontsize',16)
xticklabels({'Iso','Ortho'})
xticks([1 2])
ylabel('Fraction Suppressed')
title('First bin')

end

