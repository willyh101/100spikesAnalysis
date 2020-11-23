function plotPopRespByNumSpikes2(outVars,opts)

figure(140);clf
ensemblesToUse = opts.ensemblesToUseSpikePlot; %outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
numSpikesEachEns = outVars.numSpikesEachEns;
hzEachEns = outVars.hzEachEns; 

plotMeans = opts.plotMeansOnly;

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
    
    cmap = colorMapPicker(numel(uSpikes),outVars.defaultColorMap);
if plotMeans
   e = errorbar(uSpikes,cellfun(@(x) mean(x),data),cellfun(@(x) sem(x)*2,data));
   e.LineWidth =2;
   e.Color = rgb('black');
   hold on
   
   e2 = errorbar(0,mean(outVars.noStimPopResp(outVars.IndsUsed)),sem(outVars.noStimPopResp(outVars.IndsUsed)));
   e2.LineWidth = 2;
   e2.Color = rgb('darkgrey');
   
      xticks([0 uSpikes])
      xlim([-5 uSpikes(end)+5])

else
    p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
    
    ax=p{3};
    set(findall(gcf(),'type','line'),'markerSize',16)
    p{2}(1).Color = rgb('darkgrey');
    p{2}(2).Color = rgb('darkgrey');
    p{2}(1).LineWidth = 1;
    p{2}(2).LineWidth = 1;
end
    hold on
    r = refline(0);
    r.Color = rgb('grey');
    r.LineStyle =':';
    r.LineWidth = 2;
    title(['Ensemble Size : ' num2str(uEns(i))])
    ylabel('\Delta Pop Response')
    xlabel('Number of Added Spikes')
end