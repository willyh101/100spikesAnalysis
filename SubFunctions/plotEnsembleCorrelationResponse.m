function plotEnsembleCorrelationResponse(outVars,figNum)
if nargin<2
    figNum = 100;
end
ensSpCo = outVars.ensSpCo;
ensAlCo = outVars.ensAlCo;
ensAmCo = outVars.ensAmCo;
ensSiCo = outVars.ensSiCo;
ensNoCo = outVars.ensNoCo;


%%
ensemblesToUse = outVars.ensemblesToUse;
figure(figNum);clf
dat = {ensSpCo(ensemblesToUse), ensAlCo(ensemblesToUse), ensAmCo(ensemblesToUse),  ensSiCo(ensemblesToUse), ensNoCo(ensemblesToUse)};
names = {'Spont' 'All' 'All (v2)' 'Signal' 'Noise'};
fancyPlotSpread(dat,names);
title('Ensemble Mean Correlations by type')
ylabel('Correlation (Rho)')

%%plot Pop Response by Correlation
f3 = figure(figNum+1);
clf(f3);

popResponseEns = outVars.popResponseEns; 
numCellsEachEns = outVars.numCellsEachEns;

for i=1:5
    subplot(5,1,i)
    dataToUse = dat{i};
    scatter(dataToUse,popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
    % scatter(1:sum(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
    
    title([names{i} ' Correlation'])
    
    xlabel(['Correlation of Ensemble'])
    ylabel('Population Mean Response')
    % title('OSIs by Ensemble Size')
    set(gcf(),'Name','OSIs by Ensemble Size')
    cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
    cb.Label.String = 'Number of Cells in Ensemble';
    r = refline(0);
    r.LineStyle =':';
end

%% Plot Regression Lines

for i=1:5
    subplot(5,1,i)
    dataToUse = dat{i};
    
    dataToUse\popResponseEns(ensemblesToUse)
end
