%%
% Reproduces iso vs. ortho plot
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Should only included tuned and visually responsive cells
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig5(cellTable,cellCond)

totalNumEns = cellTable.ensNum(end);
distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

% Conditions for this analysis
cellSelectorEns = cellTable.cellEnsOSI>0.7 & cellTable.cellMeanEnsOSI>0.5; % co-tuned ensembles
cellSelectorOriIso = cellTable.cellOrisDiff == 0;
cellSelectorOriOrtho = cellTable.cellOrisDiff == 90;

cellDistDataAveOrtho=zeros(length(distBins)-1,totalNumEns);
cellDistDataAveIso=zeros(length(distBins)-1,totalNumEns);
numEnsPref = zeros(length(distBins)-1,1);
numEnsOrtho = zeros(length(distBins)-1,1);
% Loop over all ensembles
for ii = 1:totalNumEns
    % Loop over all distances
    for ll = 1:length(distBins)-1
        cellSelectorDist = cellTable.ensNum == ii & ...
            cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
        cellSelectorIso= cellSelectorDist & cellSelectorEns & cellSelectorOriIso & cellCond;
        cellSelectorOrtho= cellSelectorDist & cellSelectorEns & cellSelectorOriOrtho & cellCond;
       
        cellDistDataAveIso(ll,ii) = nanmean(cellTable.dff(cellSelectorIso));
        cellDistDataAveOrtho(ll,ii) = nanmean(cellTable.dff(cellSelectorOrtho));
        
        % Keep track of the number of ensembles used at each distance
        numEnsPref(ll) = numEnsPref(ll) + sign(sum(cellSelectorIso));
        numEnsOrtho(ll) = numEnsOrtho(ll) + sign(sum(cellSelectorOrtho));
    end
end
% Average across ensembles
respAveIso = nanmean(cellDistDataAveIso,2);
respStdErrIso = nanstd(cellDistDataAveIso,[],2)./sqrt(numEnsPref);
respAveOrtho = nanmean(cellDistDataAveOrtho,2);
respStdErrOrtho = nanstd(cellDistDataAveOrtho,[],2)./sqrt(numEnsOrtho);

% Plot the result
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


%%

% oriVals = [NaN 0:45:315];
% cellOris = oriVals(cellTable.cellPO)';
% 
% cellOrisDiff = abs(cellOris-cellTable.ensPO);
% cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
% cellOrisDiff(cellOrisDiff==135)=45;
% cellOrisDiff(cellOrisDiff==180)=0;

