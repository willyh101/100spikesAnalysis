function plotAllEnsResponse(outVars)
ensemblesToUse  = outVars.ensemblesToUse;
popResponseEns  = outVars.popResponseEns;
numCellsEachEns = outVars.numCellsEachEns;

%% Plot
f3 = figure(3);
clf(3)


%scatter(meanOSI(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
scatter(1:sum(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')

% p = polyfit(ensOSI(ensemblesToUse),popResponseEns(ensemblesToUse),1);
% f = polyval(p, ensOSI(ensemblesToUse));
% hold on
% plot(ensOSI(ensemblesToUse), f)
% hold off
xlabel('Order of being done')
ylabel('Population Mean Response')
title('Mean Response by Ensemble Size')
set(gcf(),'Name','Mean Response by Ensemble Size')
cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';