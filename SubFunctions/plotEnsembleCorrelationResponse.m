function plotEnsembleCorrelationResponse(outVars,figNum,plotRegression)
if nargin<2 || isempty(figNum)
    figNum = 100;
end
if nargin<3
    plotRegression=0;
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
    colormap(outVars.defaultColorMap)
end

%% Plot Regression Lines
if plotRegression
ensTypes = unique(numCellsEachEns(ensemblesToUse)); 

for i=1:5
    subplot(5,1,i)
    dataToUse = dat{i};
    nanFinder = isnan(dataToUse);
    dataToUse(nanFinder)=[];
    popResponseTouse = popResponseEns(ensemblesToUse);
    popResponseTouse(nanFinder)=[];
    numCellsToUse = numCellsEachEns(ensemblesToUse);
    numCellsToUse(nanFinder)=[];
    
    A = dataToUse';
    B = popResponseTouse;
    [p Rsq] = simplifiedLinearRegression(A,B);
    
    hold on
    pl = plot(dataToUse',p(1)*dataToUse+p(2),'color',rgb('dimgrey'));
    for k=1:numel(ensTypes)
        A = dataToUse(numCellsToUse==ensTypes(k));
        B = popResponseTouse(numCellsToUse==ensTypes(k))';
        [p Rsq(k+1)] = simplifiedLinearRegression(A,B);
        [x colorToUse] = colorMapPicker(numel(ensTypes),outVars.defaultColorMap,k);
        pl2(k) = plot(dataToUse',p(1)*dataToUse+p(2),'color',colorToUse);
    end

    legend([pl pl2],cellfun(@(x) num2str(x),num2cell(Rsq),'uniformoutput',0));
%     legend(pl2,num2str(Rsq,2))
end
end
