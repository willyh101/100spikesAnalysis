function plotCompareAllCellsTuning(outVars, opts)

by = opts.ensXaxis;

ensemblesToUse = outVars.ensemblesToUse;

popRespAll = outVars.allEnsResp(ensemblesToUse);
popRespIso = outVars.isoEnsResp(ensemblesToUse);
popRespOrtho = outVars.orthoEnsResp(ensemblesToUse);
popRespNotCoTuned = outVars.notCoTunedResp(ensemblesToUse);
% popRespNotTuned = outVars.notTunedEnsResp(ensemblesToUse);
popRespNotVis = outVars.notVisEnsResp(ensemblesToUse);

switch by
    case 'order'
        x = 1:numel(popRespCoTuned);
        xname = 'Order of Being Done';
    case 'osi'
        x = outVars.ensOSI(ensemblesToUse);
        xname = 'Ensemble OSI';
    case 'dist'
        x = outVars.ensGeoD(ensemblesToUse);
        xname = 'Ensemble Mean Distance';
    case 'corr'
        x = outVars.ensAlCo(ensemblesToUse);
        xname = 'Ensemble All Corr';
    case 'size'
        x = outVars.numCellsEachEns(ensemblesToUse);
        xname = 'Ensemble Size';
end

f1 = figure(200);
clf
colormap(f1, 'viridis')


subplot (1,3,1)
s1 = scatter(x, popRespIso, 'filled');
hold on
s2 = scatter(x, popRespOrtho, 'filled');
s3 = scatter(x, popRespNotCoTuned, 'filled');
% s4 = scatter(x, popRespNotTuned, 'filled');
s5 = scatter(x, popRespNotVis, 'filled');

title({'Mean Population Response','Pyramidal Cells'})
xlabel(xname)
ylabel('Population Mean Response')
% title('OSIs by Ensemble Size')
% set(gcf(),'Name','OSIs by Ensemble Size')
% % cb = colorbar('Ticks', unique(numCellsEachEns));
% % cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';
% legend([s1, s2, s3, s4, s5], {'Iso', 'Ortho', 'Not Co-Tuned', 'Not Tuned', 'Not Vis Resp'})
legend([s1, s2, s3, s5], {'Iso', 'Ortho', 'Not Co-Tuned', 'Not Vis Resp'})


subplot(1,3,2)
% 
% cats = categorical({'Iso', 'Ortho', 'Not Co-Tuned', 'Not Tuned', 'Not Vis Resp', 'All'});
% cats = reordercats(cats, {'Iso', 'Ortho', 'Not Co-Tuned', 'Not Tuned', 'Not Vis Resp', 'All'});
% data = [nanmean(popRespIso) nanmean(popRespOrtho) nanmean(popRespNotCoTuned) ...
%         nanmean(popRespNotTuned) nanmean(popRespNotVis) nanmean(popRespAll)];
% sems = [sem2(popRespIso, 2) sem2(popRespOrtho, 2) sem2(popRespNotCoTuned, 2) ...
%         sem2(popRespNotTuned, 2) sem2(popRespNotVis, 2) sem2(popRespAll, 2)];
    
cats = categorical({'Iso', 'Ortho', 'Not Co-Tuned', 'Not Vis Resp', 'All'});
cats = reordercats(cats, {'Iso', 'Ortho', 'Not Co-Tuned', 'Not Vis Resp', 'All'});
data = [nanmean(popRespIso) nanmean(popRespOrtho) nanmean(popRespNotCoTuned) ...
       nanmean(popRespNotVis) nanmean(popRespAll)];
sems = [sem2(popRespIso, 2) sem2(popRespOrtho, 2) sem2(popRespNotCoTuned, 2) ...
        sem2(popRespNotVis, 2) sem2(popRespAll, 2)];

bar(cats, data);
nbars = 1:numel(data);

hold on

er = errorbar(nbars, data, sems);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;

title({'Pop Response To Ensemble'})
ylabel('Mean Population Response')
xtickangle(45)

subplot(1,3,3)

% fancyPlotSpread({popRespIso; popRespOrtho; popRespNotCoTuned; popRespNotTuned; popRespNotVis; popRespAll}', ...
%                 {'Iso', 'Ortho', 'Not Co-Tuned', 'Not Tuned', 'Not Vis Resp', 'All'});
            
fancyPlotSpread({popRespIso; popRespOrtho; popRespNotCoTuned; popRespNotVis; popRespAll}', ...
    {'Iso', 'Ortho', 'Not Co-Tuned', 'Not Vis Resp', 'All'});
title({'Pop Response To Ensemble'})
ylabel('Mean Population Response')
xtickangle(45)

