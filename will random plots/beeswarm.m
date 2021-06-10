% function h = beeSwarmPlot(data, names)
% % for now, data should be a cell array
% 
% assert(iscell(data), 'For now, data must be a cell array by category.')
% assert(iscell(names), 'Names should be a cell array of strings.')

clear xData
ncats = numel(data);
figure(1111)
clf

sdev = 0.1;
for i=1:ncats
    hold on
    % for each category make some jittered x data
    xData{i} = sdev.*randn(numel(data{i}),1)+i;
    s = scatter(xData{i},data{i},'filled');
    s.MarkerEdgeColor = 'none';
    s.MarkerFaceAlpha = 0.6;
    hold on
    
    % then do an errorbar
    m(i) = nanmean(data{i});
    err(i) = sem(data{i});
    
end

e = errorbar(1:ncats,m,err);
e.LineWidth = 1.5;
e.Marker = 'o';
e.MarkerSize = 4;
e.MarkerFaceColor = 'k';
e.MarkerEdgeColor = 'none';
e.LineStyle = 'none';
e.Color = 'k';

xlim([0,ncats+1])


