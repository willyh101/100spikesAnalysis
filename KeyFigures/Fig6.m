function [] = Fig6(cellTable,visCells)

distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;
ensThreshs = [-inf 0.3; 0.7 inf];
meanEnsThreshs = [-inf 0.5; 0.5 inf];
spatialThresh = [-inf 400; 500 inf];

numEnsUsedTight = zeros(length(distBins)-1,2);
numEnsUsedLoose = zeros(length(distBins)-1,2);

respAveTight = zeros(length(distBins)-1,2);
respStdErrTight = zeros(length(distBins)-1,2);
respAveLoose = zeros(length(distBins)-1,2);
respStdErrLoose = zeros(length(distBins)-1,2);


% Loop over conditions
for jj = 1:2
    cellDataAveTight=[];
    cellDataAveLoose=[];
    % Loop over ensembles
    for ii = 1:160
        % Loop over distances
        for ll = 1:length(distBins)-1
            cellSelectorConds = cellTable.ensNum == ii & ...
                cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1) ...
                & cellTable.offTarget==0;
            ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(jj,1) & cellTable.cellEnsOSI<ensThreshs(jj,2)...
                & cellTable.cellMeanEnsOSI>meanEnsThreshs(jj,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(jj,2);
            
            ensSelectorSpreadTight = cellTable.cellEnsMaxD>spatialThresh(1,1) & cellTable.cellEnsMaxD<spatialThresh(1,2);
            ensSelectorSpreadLoose = cellTable.cellEnsMaxD>spatialThresh(2,1) & cellTable.cellEnsMaxD<spatialThresh(2,2);
            
            if visCells == 1
                cellSelectorTuned = cellTable.visP<0.05 & cellTable.cellOSI>0.25;
            else
                cellSelectorTuned = ones(length(cellTable.visP),1);
            end
            
            cellSelectorTight = cellSelectorConds & ensSelectorTuning & ensSelectorSpreadTight & cellSelectorTuned;
            cellSelectorLoose = cellSelectorConds & ensSelectorTuning & ensSelectorSpreadLoose & cellSelectorTuned;
            
            cellDataAveTight(ll,ii) = nanmean(cellTable.dff(cellSelectorTight));
            cellDataAveLoose(ll,ii) = nanmean(cellTable.dff(cellSelectorLoose));
            
            if sum(cellSelectorTight)~=0
                numEnsUsedTight(ll,jj) = numEnsUsedTight(ll,jj)+1;
            end
            
            if sum(cellSelectorLoose)~=0
                numEnsUsedLoose(ll,jj) = numEnsUsedLoose(ll,jj)+1;
            end
        end
    end
    
    respAveTight(:,jj) = nanmean(cellDataAveTight,2);
    respStdErrTight(:,jj) = nanstd(cellDataAveTight,[],2)./sqrt(numEnsUsedTight(:,jj));
    
    respAveLoose(:,jj) = nanmean(cellDataAveLoose,2);
    respStdErrLoose(:,jj) = nanstd(cellDataAveLoose,[],2)./sqrt(numEnsUsedLoose(:,jj));
end

colorScheme =[];
colorScheme(1,1,:) = [0 0 0]; colorScheme(1,2,:) = [0 0 0]+0.5;
colorScheme(2,1,:) = [0.494 0.184 0.556]; colorScheme(2,2,:) = [1 0 1];

figure(); clf; 
for jj = 1:2
    subplot(2,2,1+(jj-1)*2); hold on;
    plot(plotDist,respAveTight(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,1,:))
    errorbar(plotDist,respAveTight(:,jj),respStdErrTight(:,jj),'linewidth',1.5,'color',colorScheme(jj,1,:))
    plot([0 250],0*[0 250],'k--')
    set(gca,'fontsize',16)
    xlim([0 150])
    ylim([-0.1 0.1])
    
    subplot(2,2,2+(jj-1)*2); hold on;
    plot(plotDist,respAveLoose(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,2,:))
    errorbar(plotDist,respAveLoose(:,jj),respStdErrLoose(:,jj),'linewidth',1.5,'color',colorScheme(jj,2,:))
    plot([0 250],0*[0 250],'k--')
    set(gca,'fontsize',16)
    xlim([0 150])
    ylim([-0.1 0.1])
end

end

