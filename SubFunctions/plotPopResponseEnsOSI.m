function plotPopResponseEnsOSI(outVars, opts)

method = opts.ensOSImethod;

ensemblesToUse = outVars.ensemblesToUse;
popResponseEns = outVars.popResponseEns;
numCellsEachEns = outVars.numCellsEachEns;
ensOSI = outVars.(method);

figure(44);
clf

scatter(ensOSI(ensemblesToUse), popResponseEns(ensemblesToUse), [], numCellsEachEns(ensemblesToUse),'filled')
hold on

% p = polyfit(ensOSI(ensemblesToUse)',popResponseEns(ensemblesToUse),1);
% f = polyval(p, ensOSI(ensemblesToUse));
% [fs gs] = fit(ensOSI(ensemblesToUse),popResponseEns(ensemblesToUse),'poly1');
% 
% fline = plot(fs);
% fline.LineWidth = 1;
% legend('Ensemble Mean',['RSquared: ' num2str(gs.rsquare)]);


xlabel('Ensemble OSI')
ylabel('Population Mean Response')
title('OSIs by Ensemble Size')
set(gcf(),'Name','OSIs by Ensemble Size')
cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';