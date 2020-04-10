function plotEnsembleDistanceResponse(outVars,figNum,plotRegression)
if nargin<2 || isempty(figNum)
    figNum = 100;
end
if nargin<3
    plotRegression=0;
end
ensMeaD = outVars.ensMeaD;
ensGeoD = outVars.ensGeoD;
ensMaxD = outVars.ensMaxD;
ensMinD = outVars.ensMinD;

numPanels = 4;

%%
ensemblesToUse = outVars.ensemblesToUse;
figure(figNum);clf
dat = {ensMeaD(ensemblesToUse), ensGeoD(ensemblesToUse),...
    ensMaxD(ensemblesToUse),  ensMinD(ensemblesToUse)};
names = {'Mean' 'GeoMean' 'Max' 'Min'};
fancyPlotSpread(dat,names);
title('Ensemble Distance Correlations by type')
ylabel('Distance (\mum)')

%%plot Pop Response by Correlation
f4 = figure(figNum+1);
clf(f4);

popResponseEns = outVars.popResponseEns;
numCellsEachEns = outVars.numCellsEachEns;

for i=1:numPanels
    subplot(numPanels,1,i)
    dataToUse = dat{i};
    scatter(dataToUse,popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
    % scatter(1:sum(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
    hold on
    title([names{i} ' Distance'])
    
    xlabel(['Spread of Ensemble (\mum)'])
    ylabel('Population Mean Response')
    % title('OSIs by Ensemble Size')
    cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
    cb.Label.String = 'Number of Cells in Ensemble';
    r = refline(0);
    r.LineStyle =':';
    colormap(outVars.defaultColorMap);
end

%% Plot Regression Lines
if plotRegression
    ensTypes = unique(numCellsEachEns(ensemblesToUse));
    figure(f4)
    for i=1:numPanels
        subplot(numPanels,1,i)
        hold on
        
        dataToUse = dat{i};
        nanFinder = isnan(dataToUse);
        dataToUse(nanFinder)=[];
        popResponseTouse = popResponseEns(ensemblesToUse);
        popResponseTouse(nanFinder)=[];
        numCellsToUse = numCellsEachEns(ensemblesToUse);
        numCellsToUse(nanFinder)=[];
        
        A = dataToUse';
        B = popResponseTouse;
        [p Rsq pVal] = simplifiedLinearRegression(A,B);
        pV=double(pVal(1));
        pl = plot(dataToUse',p(1)*dataToUse+p(2),'color',rgb('dimgrey'));
        for k=1:numel(ensTypes)
            A = dataToUse(numCellsToUse==ensTypes(k))';
            B = popResponseTouse(numCellsToUse==ensTypes(k));
            [p Rsq(k+1), pVal] = simplifiedLinearRegression(A,B);
            pV(k+1) = double(pVal(1));
            [x colorToUse] = colorMapPicker(numel(ensTypes),outVars.defaultColorMap,k);
            pl2(k) = plot(dataToUse',p(1)*dataToUse+p(2),'color',colorToUse);
        end
        
        legend([pl pl2],cellfun(@(x) num2str(x),num2cell(pV),'uniformoutput',0));
        %     legend(pl2,num2str(Rsq,2))
    end
end
