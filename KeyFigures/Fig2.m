%%
% Reproduces the minimal distance to the ensemble plot
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig2(cellTable,cellCond)

totalNumEns = cellTable.ensNum(end);

distBins = [15:15:250];
plotDist = distBins(1:end-1) + diff(distBins(1:2))/2;

cellDistDataAve=zeros(length(distBins)-1,totalNumEns);
% Loop over all ensembles
for ii = 1:totalNumEns
    % Loop over all distances
    for ll = 1:length(distBins)-1
        cellSelectorDist = cellTable.ensNum == ii & ...
            cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
        
        cellSelector = cellSelectorDist & cellCond;
        
        cellDistDataAve(ll,ii) = nanmean(cellTable.('dff')(cellSelector));
    end
end
% Average across ensembles
respAve = nanmean(cellDistDataAve,2);
respStdErr = nanstd(cellDistDataAve,[],2)/sqrt(size(cellDistDataAve,2));

%% Plot the results
figure(); clf; hold on;
leg(1) = plot(plotDist,respAve,'k.','markersize',20);
e = errorbar(plotDist,respAve,respStdErr,'k','linewidth',2,'CapSize',0);
e.LineStyle = 'none';
plot([0 250],0*[0 250],'k--')
xlim([0 250])
xticks([0:25:250])
maxVal = abs(max(respAve+respStdErr));
minVal = abs(min(respAve-respStdErr));
ylim([-round(max(minVal,maxVal),2) round(max(minVal,maxVal),2)])
set(gca,'fontsize',16)
ylabel('Mean evoked \DeltaF/F')
xlabel('Min Dist to Stim Loc')
ylim([-0.02 0.07])

%% Fit the experimental data

expFn = 'A*exp(-x^2./sigma1)+B*exp(-x^2./sigma2)';
fnFit = fit(plotDist',respAve,expFn,'StartPoint',[0.15 -0.02 500 2e4]);
xPlot = [20:0.1:250];
leg(2) = plot(xPlot,fnFit.A*exp(-xPlot.^2./fnFit.sigma1)+fnFit.B*exp(-xPlot.^2./fnFit.sigma2),'k','linewidth',1.5);
legend(leg,{'Data','Fit'})

% sigma1 and sigma2 are really the variances of the Gaussians
fprintf('Fitted params: \n')
fprintf('A = %f, sigma1 =%f\n', fnFit.A, sqrt(fnFit.sigma1))
fprintf('B = %f, sigma2 =%f\n', fnFit.B, sqrt(fnFit.sigma2))

end

