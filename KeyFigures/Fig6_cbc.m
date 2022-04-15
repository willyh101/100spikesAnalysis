%%
% Reproduces close vs. far co-tuned plots by not averaging across all ensembles first
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function [cellResponsesTight, cellResponsesLoose] = Fig6_cbc(cellTable,cellCond)

distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

% Ensemble thresholds
ensThreshs = [-inf 0.3; 0.7 inf];
meanEnsThreshs = [-inf 0.5; 0.5 inf];
spatialThresh = [-inf 400; 500 inf];

% Conditions for this analysis
ensSelectorSpreadTight = cellTable.cellEnsMaxD>spatialThresh(1,1) & cellTable.cellEnsMaxD<spatialThresh(1,2);
ensSelectorSpreadLoose = cellTable.cellEnsMaxD>spatialThresh(2,1) & cellTable.cellEnsMaxD<spatialThresh(2,2);

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

% Plot the results
colorScheme =[];
colorScheme(1,1,:) = [0 0 0]; colorScheme(1,2,:) = [0 0 0]+0.5;
colorScheme(2,1,:) = [0.494 0.184 0.556]; colorScheme(2,2,:) = [1 0 1];
figure();
for jj = 1:2
    subplot(2,2,1+(jj-1)*2); hold on;
    plot(plotDist,respAveTight(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,1,:))
    errorbar(plotDist,respAveTight(:,jj),respStdErrTight(:,jj),'linewidth',1.5,'color',colorScheme(jj,1,:))
    plot(plotDist,respMedTight(:,jj),'*','markersize',15,'color',colorScheme(jj,1,:))
    plot([0 250],0*[0 250],'k--')
    set(gca,'fontsize',16)
    xlim([0 150])
    ylim([-0.1 0.1])
    
    subplot(2,2,2+(jj-1)*2); hold on;
    plot(plotDist,respAveLoose(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,2,:))
    errorbar(plotDist,respAveLoose(:,jj),respStdErrLoose(:,jj),'linewidth',1.5,'color',colorScheme(jj,2,:))
    plot(plotDist,respMedLoose(:,jj),'*','markersize',15,'color',colorScheme(jj,2,:))
    plot([0 250],0*[0 250],'k--')
    set(gca,'fontsize',16)
    xlim([0 150])
    ylim([-0.1 0.1])
end

%% First bin analyis

figure(); clf;
subplot(1,2,1);
boxplot([cellResponsesTight{1}; cellResponsesLoose{1}; cellResponsesTight{2}; cellResponsesLoose{2}],...
    [ones(length(cellResponsesTight{1}),1); 2*ones(length(cellResponsesLoose{1}),1); ...
    3*ones(length(cellResponsesTight{2}),1); 4*ones(length(cellResponsesLoose{2}),1)])
set(gca,'fontsize',16)
title('First bin distributions')
ylabel('\DeltaF/F')

subplot(1,2,2); hold on;
bar([sum(cellResponsesTight{1}<0)/length(cellResponsesTight{1}),...
    sum(cellResponsesLoose{1}<0)/length(cellResponsesLoose{1}),...
    sum(cellResponsesTight{2}<0)/length(cellResponsesTight{2}),...
    sum(cellResponsesLoose{2}<0)/length(cellResponsesLoose{2})])
plot([0 5],0.5+0*[0 5],'k--','linewidth',1.5)
set(gca,'fontsize',16)
title('First bin')
ylabel('Fraction Suppressed')

end