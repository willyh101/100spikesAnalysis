% function plotPopRespByStimRate(outVars)

figure(139);clf
ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
hzEachEns = outVars.hzEachEns; 

uEns = unique(numCellsEachEns(ensemblesToUse)); 
popResponseEns = outVars.popResponseEns;

clear mdat err
for i=1:numel(uEns)
    subplot(1,numel(uEns),i)
    
    rates = hzEachEns(ensemblesToUse & numCellsEachEns==uEns(i));
    uRates = unique(rates);
    
    clear data names
    for k=1:numel(uRates)
        data{k} = popResponseEns(ensemblesToUse & numCellsEachEns==uEns(i) & hzEachEns==uRates(k));
        names{k} = string(uRates(k));
    end
    mdat(i) = mean(popResponseEns(ensemblesToUse & numCellsEachEns==uEns(i) & hzEachEns==uRates(k)));
    err(i) = sem(popResponseEns(ensemblesToUse & numCellsEachEns==uEns(i) & hzEachEns==uRates(k)));
    a{i} = popResponseEns(ensemblesToUse & numCellsEachEns==uEns(i) & hzEachEns==uRates(k));
    
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

figure(299)
clf

% uRates = unique(hzEachEns(ensemblesToUse));
% clear data err
% for i=1:numel(uRates)
%     data(i) = nanmean(popResponseEns(ensemblesToUse & hzEachEns==uRates(i)));
%     err(i) = sem(popResponseEns(ensemblesToUse & hzEachEns==uRates(i)));
% end

errorbar(mdat,err)
    