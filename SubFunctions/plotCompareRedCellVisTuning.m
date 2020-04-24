function plotCompareRedCellVisTuning(outVars, opts)

by = opts.redCellXaxis;

ensemblesToUse = outVars.ensemblesToUse;
% ensOSI = outVars.ensOSI(ensemblesToUse);
popRespCoTuned = outVars.coTunedRedEnsResp(ensemblesToUse);
popRespNotCoTuned = outVars.notCoTunedRedEnsResp(ensemblesToUse);
% numCellsEachEns = outVars.numCellsEachEns(ensemblesToUse);

switch by
    case 'order'
        x = 1:numel(popRespCoTuned);
        xname = 'Order of Being Done';
    case 'osi'
        x = outVars.ensOSI(ensemblesToUse);
        xname = 'Ensemble OSI';
    case 'dist'
        x = outVars.ensMeaD(ensemblesToUse);
        xname = 'Ensemble Mean Distance';
    case 'corr'
        x = outVars.ensAlCo(ensemblesToUse);
        xname = 'Ensemble All Corr';
end

f = figure(50);
clf
colormap(f, 'viridis')

s1 = scatter(x, popRespCoTuned, 'filled');

hold on

s2 = scatter(x, popRespNotCoTuned, 'filled');
% s2.MarkerEdgeColor = 'k';

xlabel(xname)
ylabel('Red Cell Population Mean Response')
% title('OSIs by Ensemble Size')
% set(gcf(),'Name','OSIs by Ensemble Size')
% % cb = colorbar('Ticks', unique(numCellsEachEns));
% % cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';
legend([s1, s2], {'Co-Tuned', sprintf('Not Co-Tuned')})

f2 = figure(49);
clf

cats = categorical({'Co-Tuned', 'Not Co-Tuned'});
cats = reordercats(cats, {'Co-Tuned', 'Not Co-Tuned'});
data = [nanmean(popRespCoTuned) nanmean(popRespNotCoTuned)];
sems = [sem2(popRespCoTuned, 2) sem2(popRespNotCoTuned, 2)];

b = bar(cats, data);
% set(gca(), 'xticklabel', cats)

hold on

nbars = 1:numel(data);

er = errorbar(nbars, data, sems);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;

title('Pop Response To Ensemble, Red Cells')
ylabel('Mean Population Response')