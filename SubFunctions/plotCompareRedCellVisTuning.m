function plotCompareRedCellVisTuning(outVars, opts)

by = opts.redCellXaxis;

ensemblesToUse = outVars.ensemblesToUse;
% ensOSI = outVars.ensOSI(ensemblesToUse);
popRespCoTuned = outVars.coTunedRedEnsResp(ensemblesToUse);
popRespNotCoTuned = outVars.notCoTunedRedEnsResp(ensemblesToUse);
popRespNotVisResp = outVars.notVisRespRedEnsResp(ensemblesToUse);
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

f1 = figure(50);
clf
colormap(f1, 'viridis')

subplot (1,3,1)
s1 = scatter(x, popRespCoTuned, 'filled');
hold on
s2 = scatter(x, popRespNotCoTuned, 'filled');
% s2.MarkerEdgeColor = 'k';
s3 = scatter(x, popRespNotVisResp, 'filled');

title('Mean Population Response of Red Cells')
xlabel(xname)
ylabel('Red Cell Population Mean Response')
% title('OSIs by Ensemble Size')
% set(gcf(),'Name','OSIs by Ensemble Size')
% % cb = colorbar('Ticks', unique(numCellsEachEns));
% % cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';
legend([s1, s2, s3], {'Co-Tuned', 'Not Co-Tuned', sprintf('Not Vis\nResponsive')})


subplot(1,3,2)

cats = categorical({'Co-Tuned', 'Not Co-Tuned', 'Not Vis Resp'});
cats = reordercats(cats, {'Co-Tuned', 'Not Co-Tuned', 'Not Vis Resp'});
data = [nanmean(popRespCoTuned) nanmean(popRespNotCoTuned) nanmean(popRespNotVisResp)];
sems = [sem2(popRespCoTuned, 2) sem2(popRespNotCoTuned, 2) sem2(popRespNotVisResp, 2)];

bar(cats, data);
nbars = 1:numel(data);

hold on

er = errorbar(nbars, data, sems);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;

title('Pop Response To Ensemble, Red Cells')
ylabel('Mean Population Response')

subplot(1,3,3)

fancyPlotSpread([popRespCoTuned; popRespNotCoTuned; popRespNotVisResp]', {'Co-Tuned', 'Not Co-Tuned', 'Not Vis Resp'});
title('Pop Response To Ensemble, Red Cells')
ylabel('Mean Population Response')


%%----Not Red Cells----%%


ensemblesToUse = outVars.ensemblesToUse;
% ensOSI = outVars.ensOSI(ensemblesToUse);
popRespCoTuned = outVars.coTunedOtherEnsResp(ensemblesToUse);
popRespNotCoTuned = outVars.notCoTunedOtherEnsResp(ensemblesToUse);
popRespNotVisResp = outVars.notVisRespOtherEnsResp(ensemblesToUse);
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

f2 = figure(51);
clf
colormap(f2, 'viridis')

subplot (1,3,1)
s1 = scatter(x, popRespCoTuned, 'filled');
hold on
s2 = scatter(x, popRespNotCoTuned, 'filled');
% s2.MarkerEdgeColor = 'k';
s3 = scatter(x, popRespNotVisResp, 'filled');

title('Mean Population Response of Pyramidal Cells')
xlabel(xname)
ylabel('Pyramidal Cell Population Mean Response')
% title('OSIs by Ensemble Size')
% set(gcf(),'Name','OSIs by Ensemble Size')
% % cb = colorbar('Ticks', unique(numCellsEachEns));
% % cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';
legend([s1, s2, s3], {'Co-Tuned', 'Not Co-Tuned', sprintf('Not Vis\nResponsive')})


subplot(1,3,2)

cats = categorical({'Co-Tuned', 'Not Co-Tuned', 'Not Vis Resp'});
cats = reordercats(cats, {'Co-Tuned', 'Not Co-Tuned', 'Not Vis Resp'});
data = [nanmean(popRespCoTuned) nanmean(popRespNotCoTuned) nanmean(popRespNotVisResp)];
sems = [sem2(popRespCoTuned, 2) sem2(popRespNotCoTuned, 2) sem2(popRespNotVisResp, 2)];

bar(cats, data);
nbars = 1:numel(data);

hold on

er = errorbar(nbars, data, sems);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;

title('Pop Response To Ensemble, Pyramidal Cells')
ylabel('Mean Population Response')

subplot(1,3,3)

fancyPlotSpread([popRespCoTuned; popRespNotCoTuned; popRespNotVisResp]', {'Co-Tuned', 'Not Co-Tuned', 'Not Vis Resp'});
title('Pop Response To Ensemble, Pyramidal Cells')
ylabel('Mean Population Response')