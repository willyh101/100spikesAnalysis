function plotResponseByDifferenceinAnglePref_12(outVars, opts)

mRespByOriDiff = outVars.mRespByOriDiff;
sortedOSI = outVars.sortedOSI;
diffRange = -15:30:330;

figure(138);clf
% imagesc(mRespByOriDiff);
% colormap rdbu
% caxis([-0.1 0.1])
subplot(2,2,1)
meanByOriDiff = (nanmean(mRespByOriDiff));
semByOriDiff = nanstd(mRespByOriDiff)./sqrt(sum(~isnan(mRespByOriDiff)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e1.LineWidth = 1;

% plot only good OSIs
hold on
goodOSIthreshold = opts.goodOSIthresh;
datToPlot = mRespByOriDiff(sortedOSI>=goodOSIthreshold,:);
meanByOriDiff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e2 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e2.LineWidth = 1;



ylim([-0.05 0.01]) 
xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
xticks(0:30:330)
legend({['All OSIs, n= ' num2str(size(mRespByOriDiff,1))], ...
    ['OSI > ' num2str(goodOSIthreshold) ', n=' num2str(size(datToPlot,1))]})


subplot(2,2,2)
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


datToPlot = mRespByOriDiff(sortedOSI<0.25,:);
meanByOriDiff = (nanmean(datToPlot,1));
semByOriDiff = nanstd(datToPlot,[],1)./sqrt(sum(~isnan(datToPlot),1));
e = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e.LineWidth = 1;

datToPlot = mRespByOriDiff(sortedOSI>=0.25 & sortedOSI<=0.75,:);
meanByOriDiff = (nanmean(datToPlot,1));
semByOriDiff = nanstd(datToPlot,[],1)./sqrt(sum(~isnan(datToPlot),1));
e = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e.LineWidth = 1;


datToPlot = mRespByOriDiff(sortedOSI>0.75,:);
meanByOriDiff = (nanmean(datToPlot,1));
semByOriDiff = nanstd(datToPlot,[],1)./sqrt(sum(~isnan(datToPlot),1));
e = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e.LineWidth = 1;


% ylim([-0.05 0.01])
 xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')

legend('Low OSI <0.25','Mid OSI','High OSI >0.75')


subplot(2,2,3)
ensList = find(outVars.ensemblesToUse);
ensSizes = outVars.numCellsEachEns(ensList);
uEnsSizes = unique(ensSizes);

clist = colorMapPicker(numel(uEnsSizes),'viridis');
for i=1:numel(uEnsSizes);
    
    datToPlot = mRespByOriDiff(ensSizes==uEnsSizes(i),:);
    meanByOriDiff = (nanmean(datToPlot,1));
    semByOriDiff = nanstd(datToPlot,[],1)./sqrt(sum(~isnan(datToPlot),1));
    e=errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
    e.Color = clist{i};
    e.LineWidth=2;
    legendName{i} = [num2str(uEnsSizes(i)) ' n=' num2str(sum(ensSizes==uEnsSizes(i)))];
    hold on
end
legend(legendName);
 xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
title('All Ens')

subplot(2,2,4)
ensList = find(outVars.ensemblesToUse);
ensSizes = outVars.numCellsEachEns(ensList);
uEnsSizes = unique(ensSizes);

clear legendName2
c=0;
clist = colorMapPicker(numel(uEnsSizes),'viridis');
for i=1:numel(uEnsSizes);
    condToUse = ensSizes==uEnsSizes(i) & sortedOSI>=goodOSIthreshold;
    disp([num2str(uEnsSizes(i)) ' n=' num2str(sum(condToUse))]);
    if sum(condToUse)>2
        datToPlot = mRespByOriDiff(condToUse,:);
        meanByOriDiff = (nanmean(datToPlot));
        semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
        e=errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
        e.Color = clist{i};
        e.LineWidth=2;
        c=c+1;
        legendName2{c} = [num2str(uEnsSizes(i)) ' n=' num2str(sum(condToUse))];
        hold on
    elseif sum(condToUse)>0
         datToPlot = mRespByOriDiff(condToUse,:);
        meanByOriDiff = (nanmean(datToPlot,1));
        e = plot(diffRange(1:end-1)-diffRange(1),meanByOriDiff);
        e.Color = clist{i};
        e.LineWidth=2;
        c=c+1;
        legendName2{c} = [num2str(uEnsSizes(i)) ' n=' num2str(sum(condToUse))];
        hold on
    end
end
legend(legendName2);
 xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
title('Ens over OSI threshold')

figure(145)
histogram(sortedOSI)
title('Sorted OSI')

figure(132)
clf
datToPlot = mRespByOriDiff;%(sortedOSI>=goodOSIthreshold,:);
meanByOriDiff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e1.LineWidth = 3;
e1.Color = 'k';
xlabel('\Delta Preferred Angle (Deg) from Ensemble')
ylabel('\Delta Mean Population Response')
xlim([-5 185])
% ylim([-0.035 .002])

pVals = anova1(datToPlot, [], 'off');
disp(['pval ' num2str(pVals)])
box off




