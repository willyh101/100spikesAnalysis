%%
% Reproduces iso vs. ortho plot for close vs. far co-tuned condition
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Should only included tuned and visually responsive cells
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig6IsoOrtho(cellTable,cellCond)

totalNumEns = cellTable.ensNum(end);
distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

% Ensemble thresholds
ensThreshs = [0.7 inf];
meanEnsThreshs = [0.5 inf];
spatialThresh = [-inf 400; 500 inf];

% Conditions for this analysis
ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(1,1) & cellTable.cellEnsOSI<ensThreshs(1,2)...
    & cellTable.cellMeanEnsOSI>meanEnsThreshs(1,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(1,2);
cellSelectorOriIso = cellTable.cellOrisDiff==0;
cellSelectorOriOrtho = cellTable.cellOrisDiff==90;

respAveIso = zeros(length(distBins)-1,2); respAveOrtho = zeros(length(distBins)-1,2);
respStdErrIso = zeros(length(distBins)-1,2); respStdErrOrtho = zeros(length(distBins)-1,2);
numEnsUsedIso = zeros(length(distBins)-1,2); numEnsUsedOrtho = zeros(length(distBins)-1,2);

% Loop over conditions
for jj = 1:2
    
    ensSelectorSpread = cellTable.cellEnsMaxD>spatialThresh(jj,1) & cellTable.cellEnsMaxD<spatialThresh(jj,2);
    
    cellDataAveIso=zeros(length(distBins)-1,totalNumEns);
    cellDataAveOrtho=zeros(length(distBins)-1,totalNumEns);
    % Loop over ensembles
    for ii = 1:totalNumEns
        % Loop over distances
        for ll = 1:length(distBins)-1
            cellSelectorDist = cellTable.ensNum == ii & ...
                cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
                       
            cellSelectorIso = cellSelectorDist & ensSelectorTuning & ensSelectorSpread & cellSelectorOriIso & cellCond;
            cellSelectorOrtho = cellSelectorDist & ensSelectorTuning & ensSelectorSpread & cellSelectorOriOrtho & cellCond;
            
            cellDataAveIso(ll,ii) = nanmean(cellTable.dff(cellSelectorIso));
            cellDataAveOrtho(ll,ii) = nanmean(cellTable.dff(cellSelectorOrtho));
            
            % Keep track of the number of ensembles used at each distance
            numEnsUsedIso(ll,jj) = numEnsUsedIso(ll,jj)+sign(sum(cellSelectorIso));
            numEnsUsedOrtho(ll,jj) = numEnsUsedOrtho(ll,jj)+sign(sum(cellSelectorOrtho));
        end
    end
    
    % Averave across ensembles
    respAveIso(:,jj) = nanmean(cellDataAveIso,2);
    respStdErrIso(:,jj) = nanstd(cellDataAveIso,[],2)./sqrt(numEnsUsedIso(:,jj));
    
    respAveOrtho(:,jj) = nanmean(cellDataAveOrtho,2);
    respStdErrOrtho(:,jj) = nanstd(cellDataAveOrtho,[],2)./sqrt(numEnsUsedOrtho(:,jj));
end

% Plot the results
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

