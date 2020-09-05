function plotPopRespByStimRate(outVars)

figure(139);clf
ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
hzEachEns = outVars.hzEachEns; 

uEns = unique(numCellsEachEns(ensemblesToUse)); 
popResponseEns = outVars.popResponseEns;

for i=1:numel(uEns)
    subplot(1,numel(uEns),i)
    
    rates = hzEachEns(ensemblesToUse & numCellsEachEns==uEns(i));
    uRates = unique(rates);
    
    clear data names
    for k=1:numel(uRates)
        data{k} = popResponseEns(ensemblesToUse & numCellsEachEns==uEns(i) & hzEachEns==uRates(k));
        names{k} = string(uRates(k));
    end
    
    cmap = colorMapPicker(numel(uRates),outVars.defaultColorMap);

    p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
    
    ax=p{3};
    set(findall(gcf(),'type','line'),'markerSize',16)
    p{2}(1).Color = rgb('darkgrey');
    p{2}(2).Color = rgb('darkgrey');
    p{2}(1).LineWidth = 1;
    p{2}(2).LineWidth = 1;
    
    title(['Ensemble Size : ' num2str(uEns(i))])
    ylabel('\Delta Pop Response')
    xlabel('Stim Rate (Hz)')
end