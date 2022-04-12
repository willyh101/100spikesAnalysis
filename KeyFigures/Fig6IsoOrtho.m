function [] = Fig6IsoOrtho(cellTable,cellOrisDiff)

distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;
ensThreshs = [0.7 inf];
meanEnsThreshs = [0.5 inf];
spatialThresh = [-inf 400; 500 inf];

numEnsUsedIso = zeros(length(distBins)-1,2);
numEnsUsedOrtho = zeros(length(distBins)-1,2);

respAveIso = zeros(length(distBins)-1,2);
respStdErrIso = zeros(length(distBins)-1,2);
respAveOrtho = zeros(length(distBins)-1,2);
respStdErrOrtho = zeros(length(distBins)-1,2);


% Loop over conditions
for jj = 1:2
    cellDataAveIso=[];
    cellDataAveOrtho=[];
    % Loop over ensembles
    for ii = 1:160
        % Loop over distances
        for ll = 1:length(distBins)-1
            cellSelectorConds = cellTable.ensNum == ii & ...
                cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1) ...
                & cellTable.offTarget==0;
            
            ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(1,1) & cellTable.cellEnsOSI<ensThreshs(1,2)...
                & cellTable.cellMeanEnsOSI>meanEnsThreshs(1,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(1,2);
            
            ensSelectorSpread = cellTable.cellEnsMaxD>spatialThresh(jj,1) & cellTable.cellEnsMaxD<spatialThresh(jj,2);
            
            cellSelectorIso = cellOrisDiff==0 & cellTable.visP < 0.05 & cellTable.cellOSI>0.25;
            cellSelectorOrtho = cellOrisDiff==90 & cellTable.visP < 0.05 & cellTable.cellOSI>0.25;
            
            cellSelector1 = cellSelectorConds & ensSelectorTuning & ensSelectorSpread & cellSelectorIso;
            cellSelector2 = cellSelectorConds & ensSelectorTuning & ensSelectorSpread & cellSelectorOrtho;
            
            cellDataAveIso(ll,ii) = nanmean(cellTable.dff(cellSelector1));
            cellDataAveOrtho(ll,ii) = nanmean(cellTable.dff(cellSelector2));
            
            if sum(cellSelector1)~=0
                numEnsUsedIso(ll,jj) = numEnsUsedIso(ll,jj)+1;
            end
            
            if sum(cellSelector2)~=0
                numEnsUsedOrtho(ll,jj) = numEnsUsedOrtho(ll,jj)+1;
            end
        end
    end
    
    respAveIso(:,jj) = nanmean(cellDataAveIso,2);
    respStdErrIso(:,jj) = nanstd(cellDataAveIso,[],2)./sqrt(numEnsUsedIso(:,jj));
    
    respAveOrtho(:,jj) = nanmean(cellDataAveOrtho,2);
    respStdErrOrtho(:,jj) = nanstd(cellDataAveOrtho,[],2)./sqrt(numEnsUsedOrtho(:,jj));
end

%%
colorScheme =[];
colorScheme(1,:) = [0 0.447 0.741];
colorScheme(2,:) = [0.4940 0.184 0.556];

figure(); clf; 
for jj = 1:2
    subplot(2,1,jj); hold on;
    gLeg(1) = plot(plotDist,respAveIso(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(1,:));
    errorbar(plotDist,respAveIso(:,jj),respStdErrIso(:,jj),'linewidth',1.5,'color',colorScheme(1,:))
    
    gLeg(2) = plot(plotDist,respAveOrtho(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(2,:));
    errorbar(plotDist,respAveOrtho(:,jj),respStdErrOrtho(:,jj),'linewidth',1.5,'color',colorScheme(2,:))
    plot([0 250],0*[0 250],'k--')
    set(gca,'fontsize',16)
    xlim([0 150])
    
    if jj == 1
        title('Close')
    else
        title('Far')
    end
end
legend(gLeg,{'Iso','Ortho'})

end

