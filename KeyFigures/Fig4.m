%%
% Reproduces iso vs. ortho plot
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Should only included tuned and visually responsive cells
%
% Run expFigs2Through5_Dec22 to use this function
%%
function Fig4(cellTable,cellCondTuned,mouseNames)

totalNumEns = cellTable.ensNum(end);
distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

% Conditions for this analysis
cellSelectorEns = cellTable.cellEnsOSI>0.7 & cellTable.cellMeanEnsOSI>0.5; % co-tuned ensembles
cellSelectorOri = [cellTable.cellOrisDiff==0 cellTable.cellOrisDiff==45 ...
    cellTable.cellOrisDiff==90 cellTable.cellOrisDiff>0];

% Iso, 45, Ortho, Not Iso
num_conds = 4;
respAve = zeros(length(distBins)-1,2,num_conds); 
respStdErr = zeros(length(distBins)-1,2,num_conds); 

cellDistDataAve=zeros(length(distBins)-1,totalNumEns,num_conds);
numEnsUsed = zeros(length(distBins)-1,totalNumEns);
% Loop over all ensembles
for ii = 1:totalNumEns
    % Loop over all distances
    for ll = 1:length(distBins)-1
        cellSelectorDist = cellTable.ensNum == ii & ...
            cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
        for gg = 1:num_conds
            cellSelector = cellSelectorDist & cellSelectorEns ...
                        & cellSelectorOri(:,gg) & cellCondTuned;      
            cellDistDataAve(ll,ii,gg) = nanmean(cellTable.dff(cellSelector));
            % Keep track of the number of ensembles used at each distance
            numEnsUsed(ll,gg) = numEnsUsed(ll,gg)+sign(sum(cellSelector));
        end
    end
end
% Average across ensembles
for gg = 1:num_conds
    respAve(:,gg) = nanmean(cellDistDataAve(:,:,gg),2);
    respStdErr(:,gg) = nanstd(cellDistDataAve(:,:,gg),[],2)./sqrt(numEnsUsed(:,gg));
end

%% Plot the result
colorScheme = [[0,0,139]/255; [128,0,128]/255; [255,0,255]/255];
figure(231); clf;
for gg = 1:3
    subplot(1,3,gg); hold on;
    errorbar(plotDist,respAve(:,gg),respStdErr(:,gg),'linewidth',1.5,'color',colorScheme(gg,:),'capsize',0,'linewidth',2)
    plot([0 250],0*[0 250],'k--')
    xlim([0 150])
    set(gca,'fontsize',16)
    ylim([-0.1 0.125])
    if gg == 1
        title('Iso')
    elseif gg == 3
        title('Ortho')
    end
end

%%

isoData = cellDistDataAve(1,:,1);
isoData = isoData(~isnan(isoData));

orthoData = cellDistDataAve(1,:,3);
orthoData = orthoData(~isnan(orthoData));

pVal = ranksum(isoData,orthoData,'tail','right');

fprintf('Iso vs. Ortho first bin p-val: %f\n',pVal)

%%
usedEns = find(~isnan(cellDistDataAve(1,:,1))==1);
ensExpNum = [];
for ii=1:length(usedEns)
    tempExp = cellTable.expNum(cellTable.ensNum==usedEns(ii));
    ensExpNum(ii) = tempExp(1);
end

FOVs = unique(ensExpNum);
saveNames = [];
for ii = 1:length(FOVs)
    saveNames{ii} = mouseNames{FOVs(ii)};
end
mice = unique(saveNames);

fprintf('N = %d, FOVs = %d, Mice = %d\n',min(length(isoData),length(orthoData)),...
    length(FOVs),length(mice))
%%


%% Iso vs. Not Iso plot
% indicesToPlot = [1 4];
% figure(2322); 
% 
% for gg = 1:2
%     subplot(1,2,gg); hold on;
%     errorbar(plotDist,respAve(:,indicesToPlot(gg)),respStdErr(:,indicesToPlot(gg)),...
%         'linewidth',1.5,'color',colorScheme(gg,:),'capsize',0,'linewidth',2)
%     plot([0 250],0*[0 250],'k--')
%     xlim([0 150])
%     set(gca,'fontsize',16)
%     ylim([-0.1 0.125])
%     if gg == 1
%         title('Iso')
%     else
%         title('Not Iso')
%     end
% end

end


%%

% oriVals = [NaN 0:45:315];
% cellOris = oriVals(cellTable.cellPO)';
% 
% cellOrisDiff = abs(cellOris-cellTable.ensPO);
% cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
% cellOrisDiff(cellOrisDiff==135)=45;
% cellOrisDiff(cellOrisDiff==180)=0;

