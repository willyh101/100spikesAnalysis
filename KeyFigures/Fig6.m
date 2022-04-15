%%
% Reproduces close vs. far co-tuned plots
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig6(cellTable,cellCond)

totalNumEns = cellTable.ensNum(end);
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
numEnsUsedTight = zeros(length(distBins)-1,2); numEnsUsedLoose = zeros(length(distBins)-1,2);

% Loop over conditions
for jj = 1:2
    ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(jj,1) & cellTable.cellEnsOSI<ensThreshs(jj,2)...
        & cellTable.cellMeanEnsOSI>meanEnsThreshs(jj,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(jj,2);
    
    cellDataAveTight=zeros(length(distBins)-1,totalNumEns);
    cellDataAveLoose=zeros(length(distBins)-1,totalNumEns);
    % Loop over ensembles
    for ii = 1:totalNumEns
        % Loop over distances
        for ll = 1:length(distBins)-1
            cellSelectorDist = cellTable.ensNum == ii & ...
                cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
           
            cellSelectorTight = cellSelectorDist & ensSelectorTuning & ensSelectorSpreadTight & cellCond;
            cellSelectorLoose = cellSelectorDist & ensSelectorTuning & ensSelectorSpreadLoose & cellCond;
            
            cellDataAveTight(ll,ii) = nanmean(cellTable.dff(cellSelectorTight));
            cellDataAveLoose(ll,ii) = nanmean(cellTable.dff(cellSelectorLoose));
                        
            % Keep track of the number of ensembles used at each distance
            numEnsUsedTight(ll,jj) = numEnsUsedTight(ll,jj) + sign(sum(cellSelectorTight));
            numEnsUsedLoose(ll,jj) = numEnsUsedLoose(ll,jj) + sign(sum(cellSelectorLoose));
        end
    end
    
    % Averave across ensembles
    respAveTight(:,jj) = nanmean(cellDataAveTight,2);
    respStdErrTight(:,jj) = nanstd(cellDataAveTight,[],2)./sqrt(numEnsUsedTight(:,jj));
    
    respAveLoose(:,jj) = nanmean(cellDataAveLoose,2);
    respStdErrLoose(:,jj) = nanstd(cellDataAveLoose,[],2)./sqrt(numEnsUsedLoose(:,jj));
end

% Plot the results
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

