function [outInfo] = plotPopRespByNumSpikes(outVars)

figure(140);clf
ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
numSpikesEachEns = outVars.numSpikesEachEns;
hzEachEns = outVars.hzEachEns; 

uEns = unique(numCellsEachEns(ensemblesToUse)); 
popResponseEns = outVars.popResponseEns;

for i=1:numel(uEns)
    subplot(1,numel(uEns),i)
    
%     spikes = hzEachEns(ensemblesToUse & numCellsEachEns==uEns(i));
        spikes = numSpikesEachEns(ensemblesToUse & numCellsEachEns==uEns(i));

    uSpikes = unique(spikes);
    
    clear data names
    for k=1:numel(uSpikes)
        data{k} = popResponseEns(ensemblesToUse & numCellsEachEns==uEns(i) & numSpikesEachEns==uSpikes(k));
        names{k} = string(uSpikes(k));
    end
      outInfo{i}.data = data;
    outInfo{i}.names = names;
    
    cmap = colorMapPicker(numel(uSpikes),outVars.defaultColorMap);

    p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
    
    ax=p{3};
    set(findall(gcf(),'type','line'),'markerSize',16)
    p{2}(1).Color = rgb('darkgrey');
    p{2}(2).Color = rgb('darkgrey');
    p{2}(1).LineWidth = 1;
    p{2}(2).LineWidth = 1;
    
    title(['Ensemble Size : ' num2str(uEns(i))])
    ylabel('\Delta Pop Response')
    xlabel('num Spikes')
end