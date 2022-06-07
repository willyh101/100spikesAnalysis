%%
% Reproduces close vs. far co-tuned plots by not averaging across all ensembles first
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function [cellResponsesTight, cellResponsesLoose] = Fig6_cbc_dataFit(cellTable,cellCond)

% ensDistMetric = cellTable.cellEnsMaxD;
% spatialThresh = [-inf 400; 500 inf];
ensDistMetric = cellTable.cellEnsMeaD;
spatialThresh = [-inf 200; 200 inf];

binWidths = [5 10 15];
f2=[];
figure(42134); clf; hold on;

for gh = 2:2%length(binWidths)

distBins = [15:binWidths(gh):150];
plotDist = distBins(1:end-1) + binWidths(gh)/2;

% Ensemble thresholds
ensThreshs = [-inf 0.3; 0.7 inf];
meanEnsThreshs = [-inf 0.5; 0.5 inf];

% Conditions for this analysis
ensSelectorSpreadTight = ensDistMetric>spatialThresh(1,1) & ensDistMetric<spatialThresh(1,2);
ensSelectorSpreadLoose = ensDistMetric>spatialThresh(2,1) & ensDistMetric<spatialThresh(2,2);

respAveTight = zeros(length(distBins)-1,2); respAveLoose = zeros(length(distBins)-1,2);
respStdErrTight = zeros(length(distBins)-1,2); respStdErrLoose = zeros(length(distBins)-1,2);
respMedTight = zeros(length(distBins)-1,2); respMedLoose = zeros(length(distBins)-1,2);

% Loop over conditions
for jj = 1:2
    ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(jj,1) & cellTable.cellEnsOSI<ensThreshs(jj,2)...
        & cellTable.cellMeanEnsOSI>meanEnsThreshs(jj,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(jj,2);
    % Loop over distances
    for ll = 1:length(distBins)-1
        cellSelectorDist = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);

        cellSelectorTight = cellSelectorDist & ensSelectorTuning & ensSelectorSpreadTight & cellCond;
        cellSelectorLoose = cellSelectorDist & ensSelectorTuning & ensSelectorSpreadLoose & cellCond;
                
        respAveTight(ll,jj) = nanmean(cellTable.dff(cellSelectorTight));
        respStdErrTight(ll,jj) = nanstd(cellTable.dff(cellSelectorTight))./sqrt(sum(cellSelectorTight));
        respMedTight(ll,jj)=nanmedian(cellTable.dff(cellSelectorTight));
        
        respAveLoose(ll,jj) = nanmean(cellTable.dff(cellSelectorLoose));
        respStdErrLoose(ll,jj) = nanstd(cellTable.dff(cellSelectorLoose))./sqrt(sum(cellSelectorLoose));
        respMedLoose(ll,jj)=nanmedian(cellTable.dff(cellSelectorLoose));
        
        % Save the first bin
        if ll == 1
            cellResponsesTight{jj} = cellTable.dff(cellSelectorTight);
            cellResponsesLoose{jj} = cellTable.dff(cellSelectorLoose);
        end
    end
end

expFn = 'A*exp(-x^2./sigma1)+B*exp(-x^2./sigma2)';
if gh == 1
    f2{gh} = fit(plotDist',respAveTight(:,2),expFn,'StartPoint',[0.15 -0.02 200 2e4]);
else
    f2{gh} = fit(plotDist',respAveTight(:,2),expFn,'StartPoint',[0.15 -0.02 50 2e4]);
end
g2{gh} = fit(plotDist',respAveLoose(:,2),expFn,'StartPoint',[0.15 -0.02 500 2e4]);

xPlot = [15:0.1:plotDist(end)];
colorScheme =[];
colorScheme(1,1,:) = [0 0 0]; colorScheme(1,2,:) = [0 0 0]+0.5;
colorScheme(2,1,:) = [0.494 0.184 0.556]; colorScheme(2,2,:) = [1 0 1];
plot(plotDist,respAveTight(:,2),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(2,1,:))
plot(xPlot,f2{gh}.A*exp(-xPlot.^2./f2{gh}.sigma1)+f2{gh}.B*exp(-xPlot.^2./f2{gh}.sigma2),'linewidth',1.5,'color',colorScheme(2,1,:));

plot(plotDist,respAveLoose(:,2),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(2,2,:))
plot(xPlot,g2{gh}.A*exp(-xPlot.^2./g2{gh}.sigma1)+g2{gh}.B*exp(-xPlot.^2./g2{gh}.sigma2),'linewidth',1.5,'color',colorScheme(2,2,:));

plot([0 150],0*[0 150],'k--')
set(gca,'fontsize',16)
xlabel('Min Dist')
ylabel('dF/F')
legend({'Data Tight','Fit Tight','Data Loose','Fit Loose'}) 

% % Plot the results
% colorScheme =[];
% colorScheme(1,1,:) = [0 0 0]; colorScheme(1,2,:) = [0 0 0]+0.5;
% colorScheme(2,1,:) = [0.494 0.184 0.556]; colorScheme(2,2,:) = [1 0 1];
% figure();
% for jj = 1:2
%     subplot(2,2,1+(jj-1)*2); hold on;
%     plot(plotDist,respAveTight(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,1,:))
%     errorbar(plotDist,respAveTight(:,jj),respStdErrTight(:,jj),'linewidth',1.5,'color',colorScheme(jj,1,:))
%     plot(plotDist,respMedTight(:,jj),'*','markersize',15,'color',colorScheme(jj,1,:))
%     plot([0 250],0*[0 250],'k--')
%     set(gca,'fontsize',16)
%     xlim([0 150])
%     ylim([-0.1 0.1])
%     
%     subplot(2,2,2+(jj-1)*2); hold on;
%     plot(plotDist,respAveLoose(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,2,:))
%     errorbar(plotDist,respAveLoose(:,jj),respStdErrLoose(:,jj),'linewidth',1.5,'color',colorScheme(jj,2,:))
%     plot(plotDist,respMedLoose(:,jj),'*','markersize',15,'color',colorScheme(jj,2,:))
%     plot([0 250],0*[0 250],'k--')
%     set(gca,'fontsize',16)
%     xlim([0 150])
%     ylim([-0.1 0.1])
% end

end

%%
'here';

ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(2,1) & cellTable.cellEnsOSI<ensThreshs(2,2)...
    & cellTable.cellMeanEnsOSI>meanEnsThreshs(2,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(2,2);

cellSelectorTight =  ensSelectorTuning & ensSelectorSpreadTight & cellCond ...
    & ~isnan(cellTable.dff);
cellSelectorLoose =  ensSelectorTuning & ensSelectorSpreadLoose & cellCond ...
    & ~isnan(cellTable.dff);

expFn = 'A*exp(-x^2./sigma1)+B*exp(-x^2./sigma2)';
fTight = fit(cellTable.cellDist(cellSelectorTight),cellTable.dff(cellSelectorTight),expFn,'StartPoint',[0.15 -0.02 100 2e4]);

fLoose = fit(cellTable.cellDist(cellSelectorLoose),cellTable.dff(cellSelectorLoose),expFn,'StartPoint',[0.15 -0.02 100 2e4]);

figure(4232888); clf; hold on;
% plot(cellTable.cellDist(cellSelectorTight),cellTable.dff(cellSelectorTight),'.')
xPlot = [15:0.1:plotDist(end)];
plot(xPlot,fTight.A*exp(-xPlot.^2./fTight.sigma1)+fTight.B*exp(-xPlot.^2./fTight.sigma2),'linewidth',1.5);
plot(xPlot,fLoose.A*exp(-xPlot.^2./fLoose.sigma1)+fLoose.B*exp(-xPlot.^2./fLoose.sigma2),'linewidth',1.5);

plot([0 150],0*[0 150],'k--')

    

end