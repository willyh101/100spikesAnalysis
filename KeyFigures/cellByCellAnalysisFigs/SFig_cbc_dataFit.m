%%
% Reproduces the minimal distance to the ensemble plot by not averaging
% across all ensembles first
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function SFig_cbc_dataFit(cellTable,cellCond)


colorScheme = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880];

binWidths = [5 10 15 20];
figure(32177); clf;
f2=[];
zero_crossing=[];

for ii = 1:length(binWidths)

distBins = [15:binWidths(ii):250];
plotDist = distBins(1:end-1) + binWidths(ii)/2;

cellDistDataAve = zeros(1,length(distBins)-1);
cellDistDataMedian = zeros(1,length(distBins)-1);
cellDistDataErr = zeros(1,length(distBins)-1);
% Loop over all distances
for ll = 1:length(distBins)-1
    cellSelectorDist = cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
    cellSelector = cellSelectorDist & cellCond;
    
    % Average across distances
    cellDistDataAve(ll) = nanmean(cellTable.dff(cellSelector));
    cellDistDataMedian(ll) = nanmedian(cellTable.dff(cellSelector));
    cellDistDataErr(ll) = nanstd(cellTable.dff(cellSelector))/sqrt(sum(cellSelector));
end

subplot(1,length(binWidths)+1,ii+1); hold on;
hold on
leg(1) = plot(plotDist,cellDistDataAve,'k.-','markersize',16);
errorbar(plotDist,cellDistDataAve,cellDistDataErr,'k','linewidth',1.5)
plot([0 distBins(end)],0*[0 distBins(end)],'k--');
set(gca,'fontsize',16)
xlabel('Min Dist to Stim Location')
ylabel('Mean Evoked \DeltaF/F')
xlim([0 250])
xticks([0:50:250])

xPlot = [15:0.1:plotDist(end)];
expFn = 'A*exp(-x^2./sigma1)+B*exp(-x^2./sigma2)';
f2{ii+1} = fit(plotDist',cellDistDataAve',expFn,'StartPoint',[0.15 -0.02 500 2e4]);
leg(2) = plot(xPlot,f2{ii+1}.A*exp(-xPlot.^2./f2{ii+1}.sigma1)+f2{ii+1}.B*exp(-xPlot.^2./f2{ii+1}.sigma2),'linewidth',1.5,'color',colorScheme(ii+1,:));
plot([0 250],0*[0 250],'k--')
myfun = @(x) f2{ii+1}.A*exp(-x.^2./f2{ii+1}.sigma1)+f2{ii+1}.B*exp(-x.^2./f2{ii+1}.sigma2);
zero_crossing(ii+1) = fzero(myfun,25);
maxmin_ratio = max(myfun(xPlot))/abs(min(myfun(xPlot)));

legend(leg,{'Binned Data','Data Fit'})

fprintf('Zero-crossing: %.2f, Max/Min: %.2f\n',zero_crossing(ii+1),maxmin_ratio)

title(sprintf('Bin width: %.2f', binWidths(ii)));
end


% [maxVal,maxIndex]=max(f2{ii+1}.A*exp(-xPlot.^2./f2{ii+1}.sigma1)+f2{ii+1}.B*exp(-xPlot.^2./f2{ii+1}.sigma2));
% [minVal,minIndex]=min(f2{ii+1}.A*exp(-xPlot.^2./f2{ii+1}.sigma1)+f2{ii+1}.B*exp(-xPlot.^2./f2{ii+1}.sigma2));
% 
% plot(xPlot(maxIndex),maxVal,'s','color',colorScheme(4,:),'MarkerFaceColor',colorScheme(4,:),'markersize',10)
% plot(xPlot(minIndex),minVal,'s','color',colorScheme(4,:),'MarkerFaceColor',colorScheme(4,:),'markersize',10)
% plot(zero_crossing(4),0,'.','color',colorScheme(4,:),'markersize',20)
% xlim([0 150])
% legend(leg,{'Binned Data','Data Fit'})

%%
subplot(1,length(binWidths)+1,1); hold on;

cellCond2 = cellCond & ~isnan(cellTable.dff);
plot(cellTable.cellDist(cellCond2),cellTable.dff(cellCond2),'k.','markersize',1)
ylim([-4 5])
ylabel('Mean Evoked \DeltaF/F')

set(gca,'colororder',colorScheme(1,:))
yyaxis right

f2{1} = fit(cellTable.cellDist(cellCond2),cellTable.dff(cellCond2),expFn,'StartPoint',[0.15 -0.02 500 2e4]);
leg(2) = plot(xPlot,f2{1}.A*exp(-xPlot.^2./f2{1}.sigma1)+f2{1}.B*exp(-xPlot.^2./f2{1}.sigma2),...
    'linewidth',1.5,'color',colorScheme(1,:));
plot([0 250],0*[0 250],'k--')
myfun = @(x) f2{1}.A*exp(-x.^2./f2{1}.sigma1)+f2{1}.B*exp(-x.^2./f2{1}.sigma2);
zero_crossing(1) = fzero(myfun,25);
maxmin_ratio = max(myfun(xPlot))/abs(min(myfun(xPlot)));
set(gca,'fontsize',16)
xlabel('Min Dist to Stim Location')
xlim([0 250])
xticks([0:50:250])
ylim([-0.16 0.2])

fprintf('Zero-crossing: %.2f, Max/Min: %.2f\n',zero_crossing(1),maxmin_ratio)

% legend(leg,{'Data Fit'})

%%

figure(888);clf; hold on;

for ii = 1:4
legH(ii) = plot(xPlot,f2{ii}.A*exp(-xPlot.^2./f2{ii}.sigma1)+f2{ii}.B*exp(-xPlot.^2./f2{ii}.sigma2),...
    'linewidth',1.5,'color',colorScheme(ii,:));

plot(zero_crossing(ii),0,'.','color',colorScheme(ii,:),'markersize',20)

[maxVal,maxIndex]=max(f2{ii}.A*exp(-xPlot.^2./f2{ii}.sigma1)+f2{ii}.B*exp(-xPlot.^2./f2{ii}.sigma2));
[minVal,minIndex]=min(f2{ii}.A*exp(-xPlot.^2./f2{ii}.sigma1)+f2{ii}.B*exp(-xPlot.^2./f2{ii}.sigma2));

plot(xPlot(maxIndex),maxVal,'s','color',colorScheme(ii,:),'MarkerFaceColor',colorScheme(ii,:),'markersize',10)
plot(xPlot(minIndex),minVal,'s','color',colorScheme(ii,:),'MarkerFaceColor',colorScheme(ii,:),'markersize',10)

fprintf('Narrow: %.2f, Broad: %.2f\n',sqrt(f2{ii}.sigma1),sqrt(f2{ii}.sigma2))

end
plot([0 distBins(end)],0*[0 distBins(end)],'k--');
set(gca,'fontsize',16)
xlabel('Min Dist to Stim Location')
xlim([0 100])
xticks([0:50:250])
legend(legH,{'Unbinned data','Bin width=5','Bin width=10','Bin width=15'})
ylabel('Mean Evoked \DeltaF/F')
title('Line fits')


end

