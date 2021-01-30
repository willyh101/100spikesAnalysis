function plotResponseByDiffAnglePrefPosNeg_12(outVars,opts)

mRespByOriDiffPos = outVars.mRespByOriDiffPos;
mRespByOriDiffNeg = outVars.mRespByOriDiffNeg;
sortedOSI = outVars.sortedOSI;
diffRange = -15:30:330;

figure(140);clf
subplot(3,2,1)
meanByOriDiff = (nanmean(mRespByOriDiffPos));
semByOriDiff = nanstd(mRespByOriDiffPos)./sqrt(sum(~isnan(mRespByOriDiffPos)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e1.LineWidth = 1;

% plot only good OSIs
hold on
goodOSIthreshold = opts.goodOSIthresh;
datToPlot = mRespByOriDiffPos(sortedOSI>=goodOSIthreshold,:);
meanByOriDiff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e2 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e2.LineWidth = 1;
title('Change in Cells that are excited')
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
legend({['All OSIs, n= ' num2str(size(mRespByOriDiffPos,1))], ...
    ['OSI > ' num2str(goodOSIthreshold) ', n=' num2str(size(datToPlot,1))]})


%neg
subplot(3,2,2)
meanByOriDiff = (nanmean(mRespByOriDiffNeg));
semByOriDiff = nanstd(mRespByOriDiffNeg)./sqrt(sum(~isnan(mRespByOriDiffNeg)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e1.LineWidth = 1;

% plot only good OSIs
hold on
goodOSIthreshold = opts.goodOSIthresh;
datToPlot = mRespByOriDiffNeg(sortedOSI>=goodOSIthreshold,:);
meanByOriDiff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e2 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e2.LineWidth = 1;

title('Change in Cells that are Inhibited')
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
legend({['All OSIs, n= ' num2str(size(mRespByOriDiffNeg,1))], ...
    ['OSI > ' num2str(goodOSIthreshold) ', n=' num2str(size(datToPlot,1))]})

%% Plot by size
subplot(3,2,3)
ensList = find(outVars.ensemblesToUse);
ensSizes = outVars.numCellsEachEns(ensList);
uEnsSizes = unique(ensSizes);

clist = colorMapPicker(numel(uEnsSizes),'viridis');
for i=1:numel(uEnsSizes);
    
    datToPlot = mRespByOriDiffPos(ensSizes==uEnsSizes(i),:);
    meanByOriDiff = (nanmean(datToPlot));
    semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
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
title('All Ens Pos')

subplot(3,2,4)
for i=1:numel(uEnsSizes);
    
    datToPlot = mRespByOriDiffNeg(ensSizes==uEnsSizes(i),:);
    meanByOriDiff = (nanmean(datToPlot));
    semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
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
title('All Ens Neg')

%% and filter by good OSI
clist = colorMapPicker(numel(uEnsSizes),'viridis');
subplot(3,2,5)
clear legendName2
c=0;
clist = colorMapPicker(numel(uEnsSizes),'viridis');
for i=1:numel(uEnsSizes);
    condToUse = ensSizes==uEnsSizes(i) & sortedOSI>=goodOSIthreshold;
    disp([num2str(uEnsSizes(i)) ' n=' num2str(sum(condToUse))]);
    if sum(condToUse)>2
        datToPlot = mRespByOriDiffPos(condToUse,:);
        meanByOriDiff = (nanmean(datToPlot));
        semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
        e=errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
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
title('Ens over OSI threshold Pos')

subplot(3,2,6)
clear legendName2
c=0;
clist = colorMapPicker(numel(uEnsSizes),'viridis');
for i=1:numel(uEnsSizes);
    condToUse = ensSizes==uEnsSizes(i) & sortedOSI>=goodOSIthreshold;
    disp([num2str(uEnsSizes(i)) ' n=' num2str(sum(condToUse))]);
    if sum(condToUse)>2
        datToPlot = mRespByOriDiffNeg(condToUse,:);
        meanByOriDiff = (nanmean(datToPlot));
        semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
        e=errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
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
title('Ens over OSI thresholdNeg')


