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
colormap(f, 'viridis')


subplot (1,3,1)
s1 = scatter(x, popRespRedVR, 'filled');
hold on
s2 = scatter(x, popRespRedNVR, 'filled');
% s2.MarkerEdgeColor = 'k';

title('Mean Population Response of Red Cells')
xlabel(xname)
ylabel('Red Cell Population Mean Response')
% title('OSIs by Ensemble Size')
% set(gcf(),'Name','OSIs by Ensemble Size')
% % cb = colorbar('Ticks', unique(numCellsEachEns));
% % cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';
legend([s1, s2], {'Visually Responsive', sprintf('Not Visually\nResponsive')})


subplot(1,3,2)

cats = categorical({'Vis Resp', 'Not Vis Resp'});
cats = reordercats(cats ,{'Vis Resp', 'Not Vis Resp'});
data = [mean(popRespRedVR) mean(popRespRedNVR)];
sems = [sem2(popRespRedVR, 2) sem2(popRespRedNVR, 2)];

b = bar(cats, data);
nbars = 1:numel(data);

hold on

er = errorbar(nbars, data, sems);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;

title('Pop Response To Ensemble, Red Cells')
ylabel('Mean Population Response')


subplot(1,3,3)

f = fancyPlotSpread([popRespRedVR; popRespRedNVR]', {'Vis Resp', 'Not Vis Resp'});
title('Pop Response To Ensemble, Red Cells')
ylabel('Mean Population Response')