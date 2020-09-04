function plotRedandPyrConnTogether(outVars, opts)

mRespByOriDiffNotRed = outVars.mRespByOriDiffNotRed;
mRespByOriDiffRed = outVars.mRespByOriDiffRed;

sortedOSI = outVars.sortedOSIred;
diffRange = -22.5:45:315;

goodOSIthreshold = opts.goodOSIthresh;
redCellName = opts.redCellName;

figure(177);clf

% pyramids
datToPlot = mRespByOriDiffNotRed(sortedOSI>=goodOSIthreshold,:);
meanByOriDiff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e1.LineWidth = 2;

hold on

% reds
datToPlot = mRespByOriDiffRed(sortedOSI>=goodOSIthreshold,:);
meanByOriDiff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e2 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
e2.LineWidth = 2;


% ylim([-0.05 0.01]) 
xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Mean Population Response')
xticks(0:45:135)
legend('Pyramidal Cells', redCellName)
title('Population Response to Ensemble Stim, by tuning')

figure(178);
clf

osisToUse = [0 .1 .2 .3 .4 .5 .6 .7];

for i=1:numel(osisToUse)
    osiToPlot = osisToUse(i);
    
    subplot(2,4,i)
    
    datToPlot = mRespByOriDiffNotRed(sortedOSI>=osiToPlot,:);
    nEns = size(datToPlot,1);
    meanByOriDiff = (nanmean(datToPlot));
    semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
    e = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
    e.LineWidth = 1.5;
    
    pValNotRed = anova1(datToPlot, [], 'off');
    disp(['pval not red ' num2str(pValNotRed)])

    hold on
    % reds
    datToPlot = mRespByOriDiffRed(sortedOSI>=osiToPlot,:);
    meanByOriDiff = (nanmean(datToPlot));
    semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
    e = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDiff,semByOriDiff);
    e.LineWidth = 1.5;
    
    pValNotRed = anova1(datToPlot, [], 'off');
    disp(['pval red ' num2str(pValNotRed)])
    
    title(['OSI > ' num2str(osiToPlot) ' (n=' num2str(nEns) ' ensembles)'])
    xticks(0:45:180)
    if i ==1 || i==5
        ylabel('\Delta Mean Population Response')
    end
    if i > 4
        xlabel('\Delta Preferred Angle (Deg)')
    end
    
    

    

end
legend('Pyramidal Cells', redCellName)
    