%%
% Reproduces the minimal distance to the ensemble plot
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function SFig_dataFit(cellTable,cellCond)

totalNumEns = cellTable.ensNum(end);

binWidths = [5 10 15 20];
figure(32177); clf;
f2=[];
zero_crossing=[];

for gh = 1:length(binWidths)


distBins = [15:binWidths(gh):250];
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
%%
% Plot the results
subplot(1,4,gh); hold on;
leg(1) = plot(plotDist,respAve,'k.','markersize',20);
e = errorbar(plotDist,respAve,respStdErr,'k','linewidth',2,'CapSize',0);
e.LineStyle = 'none';
plot([0 250],0*[0 250],'k--')
xlim([0 250])
xticks([0:50:250])
maxVal = abs(max(respAve+respStdErr));
minVal = abs(min(respAve-respStdErr));
ylim([-round(max(minVal,maxVal),2) round(max(minVal,maxVal),2)])
set(gca,'fontsize',16)
ylabel('Mean evoked \DeltaF/F')
xlabel('Min Dist to Stim Loc')
ylim([-0.02 0.15])

title(sprintf('Bin width: %.2f', binWidths(gh)));

%%

expFn = 'A*exp(-x^2./sigma1)+B*exp(-x^2./sigma2)';
f2{gh} = fit(plotDist',respAve,expFn,'StartPoint',[0.15 -0.02 500 2e4]);
xPlot = [plotDist(1):0.01:250];
leg(2) = plot(xPlot,f2{gh}.A*exp(-xPlot.^2./f2{gh}.sigma1)+f2{gh}.B*exp(-xPlot.^2./f2{gh}.sigma2),'k','linewidth',1.5);

myfun = @(x) f2{gh}.A*exp(-x.^2./f2{gh}.sigma1)+f2{gh}.B*exp(-x.^2./f2{gh}.sigma2);
maxmin_ratio = max(myfun(xPlot))/abs(min(myfun(xPlot)));

zero_crossing(gh) = fzero(myfun,25);
fprintf('Zero-crossing: %.2f, Max/Min: %.2f\n',zero_crossing(gh),maxmin_ratio)

legend(leg,{'Data','Fit'})
end

end

