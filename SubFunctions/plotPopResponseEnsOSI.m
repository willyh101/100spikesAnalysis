function [results ] = plotPopResponseEnsOSI(outVars, opts)

if isfield(opts,'plotFit')
    plotFit = opts.plotFit;
else
    plotFit = 0;
end

if isfield(opts,'figToUse1')
    figToUse = opts.figToUse1;
else
    figToUse = 44;
end

method = opts.ensOSImethod;

ensemblesToUse = outVars.ensemblesToUse;
popResponseEns = outVars.popResponseEns;
numCellsEachEns = outVars.numCellsEachEns;
ensOSI = outVars.(method);

figure(figToUse);
clf

scatter(ensOSI(ensemblesToUse), popResponseEns(ensemblesToUse), [], numCellsEachEns(ensemblesToUse),'filled')
hold on

results=[];
if plotFit
%     p = polyfit(ensOSI(ensemblesToUse)',popResponseEns(ensemblesToUse),1);
%     f = polyval(p, ensOSI(ensemblesToUse));
    x = ensOSI(ensemblesToUse)';
    y = popResponseEns(ensemblesToUse);
    nanEither = isnan(x) | isnan(y); 
    
    [fs, gs] = fit(x(~nanEither),y(~nanEither),'poly1');
    
    [p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither));
    
    fline = plot(fs);
    fline.LineWidth = 1;
%     legend('Ensemble Mean',['RSquared: ' num2str(gs.rsquare)]);
    legend('Ensemble Mean',['RSquared: ' num2str(gs.rsquare) '. pVal: ' num2str(pVal(1))]);
%     
     results.p = p;
     results.Rsq = Rsq;
     results.pVal = pVal;
    results.fs = fs;
    results.gs = gs;
end

xlabel('Ensemble OSI')
ylabel('Population Mean Response')
title(['Resp by ' method])
set(gcf(),'Name','OSIs by Ensemble Size')
cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';