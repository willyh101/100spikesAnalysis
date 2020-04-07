function plotResponseBySize(outVars)
ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
popResponseEns = outVars.popResponseEns;
noStimPopResp = outVars.noStimPopResp;


%% more simple, take the means, population response by ensemble size
clear avg err ns ens2plt
f6 = figure(6);
clf(f6)
numEns = numel(unique(numCellsEachEns(ensemblesToUse)));
uniqueEns = unique(numCellsEachEns(ensemblesToUse));

x = 1:numEns;
clear data names
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{i} = popResponseEns(ens2plot);
    names{i} = string(uniqueEns(i));
    avg(i) = mean(popResponseEns(ens2plot));
    err(i) = sem(popResponseEns(ens2plot));
    ns(i) = numel(popResponseEns(ens2plot));
end

data{end+1} = noStimPopResp;
names{end+1} = 'No Stim';

cmap=colormap(viridis(numEns));
cmap(end+1,:)=rgb('grey');
p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
% bar(x, avg)
% hold on
% er = errorbar(x, avg, err);
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% hold off
ylabel('Population Response (vis responsive)')
% xticklabels(uniqueEns)
% xticks = 1:6;
title('Mean population response to holo')
xlabel('Ensemble Size')
set(gcf(),'Name','Mean population response to holo')
% ns

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

 r = refline(0);
    r.LineStyle=':';
    r.Color = rgb('grey');

pValEnselbeSize = anovan(popResponseEns(ensemblesToUse),numCellsEachEns(ensemblesToUse)','display','off');

disp(['Anova pVal between sizes: ' num2str(pValEnselbeSize)]);
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))
uniqueEns(end+1) = 0; 
for i=1:size(data,2)
    %     prs = ranksum(data{i},0);
    psr = signrank(data{i});
    
    [h p ] = ttest(data{i},0);
    disp(['Size: ' num2str(uniqueEns(i)) ', n=' num2str(numel(data{i}))...
        '. Mean: ' num2str(nanmean(data{i}),2) ' ' char(177) ' ' num2str(sem(data{i}),2)...
        '. Signed Rank: ' num2str(psr,2) '. ttest: ' num2str(p,2)])
end