function p = fancyPlotSpread(data, names)
% wrapper function for plotSpread
% sets my preferred defaults

cmap = colormap(viridis(numel(names)));
p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors', cmap);

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1.5;
p{2}(2).LineWidth = 1.5;
uistack(p{2},'bottom')
