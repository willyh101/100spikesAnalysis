function f = plotCompareRedCellVisResp(outVars, opts)

by = opts.redCellXaxis;

ensemblesToUse = outVars.ensemblesToUse;
% ensOSI = outVars.ensOSI(ensemblesToUse);
popRespRedVR = outVars.popRespEnsRedVisResp(ensemblesToUse);
popRespRedNVR = outVars.popRespEnsRedNotVisResp(ensemblesToUse);
% numCellsEachEns = outVars.numCellsEachEns(ensemblesToUse);

switch by
    case 'order'
        x = 1:numel(popRespRedVR);
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

f = figure(47);
clf
colormap(f, 'jet')

s1 = scatter(x, popRespRedVR, 'filled');

hold on

s2 = scatter(x, popRespRedNVR, 'filled');
% s2.MarkerEdgeColor = 'k';

xlabel(xname)
ylabel('Red Cell Population Mean Response')
% title('OSIs by Ensemble Size')
% set(gcf(),'Name','OSIs by Ensemble Size')
% % cb = colorbar('Ticks', unique(numCellsEachEns));
% % cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';
legend([s1, s2], {'Visually Responsive', sprintf('Not Visually\nResponsive')})

f2 = figure(48);
clf

cats = categorical({'Vis Resp', 'Not Vis Resp'});
cats = reordercats(cats ,{'Vis Resp', 'Not Vis Resp'});
data = [mean(popRespRedVR) mean(popRespRedNVR)];
sems = [sem2(popRespRedVR, 2) sem2(popRespRedNVR, 2)];

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