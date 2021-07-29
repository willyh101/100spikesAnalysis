function plotResponseByDifferenceinAnglePrefRed_12(outVars, opts)

mRespByOriDiff = outVars.mRespByOriDiffRed;
sortedOSI = outVars.sortedOSIred;
diffRange = -15:30:330;

figure(144);clf
% imagesc(mRespByOriDiff);
% colormap rdbu
% caxis([-0.1 0.1])
subplot(1,2,1)
meanByOriDiff = (nanmean(mRespByOriDiff));
semByOriDiff = nanstd(mRespByOriDiff)./sqrt(sum(~isnan(mRespByOriDiff)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e1.LineWidth = 2;

% plot only good OSIs
hold on
goodOSIthreshold = opts.goodOSIthresh;
datToPlot = mRespByOriDiff(sortedOSI>=goodOSIthreshold,:);
meanByOriDiff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e2 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e2.LineWidth = 2;



% ylim([-0.05 0.01]) 
xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
xticks(0:30:160)
legend({['All OSIs, n= ' num2str(size(mRespByOriDiff,1))], ...
    ['OSI > ' num2str(goodOSIthreshold) ', n=' num2str(size(datToPlot,1))]})


subplot(1,2,2)
% imToPlot = mRespByOriDiff;
% imToPlot(isnan(imToPlot))=0;
% imagesc(imToPlot)
% caxis([-0.2 0.2])
% xlim([0.5 7.5])
% pT1 = prctile(sortedOSI,25);
% pT2 = prctile(sortedOSI,50);
% pT3 = prctile(sortedOSI,75);

numEns = numel(sortedOSI); 
hold on
% datToPlot = mRespByOriDiff(1:round(numEns/4),:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(round(numEns/4):round(numEns/2),:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(round(numEns/2):round(numEns*3/4),:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(round(numEns*3/4):end,:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)

% 
% datToPlot = mRespByOriDiff(sortedOSI<0.25,:);
% meanByOriDiff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(sortedOSI>=0.25 & sortedOSI<=0.75,:);
% meanByOriDiff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(sortedOSI>0.75,:);
% meanByOriDiff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff)


% ylim([-0.05 0.01])
 xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')

legend('Low OSI <0.25','Mid OSI','High OSI >0.75')

figure(147)
histogram(sortedOSI)
title('Sorted OSI')

figure(1444)
clf
osiToPlot = 0.5;
meanByOriDiff = (nanmean(mRespByOriDiff(mRespByOriDiff >= osiToPlot)));
semByOriDiff = nanstd(mRespByOriDiff)./sqrt(sum(~isnan(mRespByOriDiff)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e1.LineWidth = 2;
xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
xticks(0:30:160)



