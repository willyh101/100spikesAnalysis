ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
popResponseEns = outVars.popResponseEns;

    if plotAllNoStim
noStimPopResp = outVars.noStimPopResp;
    else
       noStimPopResp =  outVars.noStimPopResp(outVars.IndsUsed);
    end

%%more simple, take the means, population response by ensemble size
clear avg err ns ens2plt
f6 = figure(6);
clf(f6)
numEns = numel(unique(numCellsEachEns(ensemblesToUse)));
uniqueEns = unique(numCellsEachEns(ensemblesToUse));

x = 1:numEns;
clear data names cmap
c=1;
for i=1:numEns
    c=c+1;
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{c} = popResponseEns(ens2plot);
    names{c} = string(uniqueEns(i));
    avg(c) = mean(popResponseEns(ens2plot));
    err(c) = sem(popResponseEns(ens2plot));
    ns(c) = numel(popResponseEns(ens2plot));
end

data{1} = noStimPopResp;
names{1} = 'Control';
names{2} = sprintf('Stimulation');

% cmap = colormap(outVars.defaultColorMap);
% cmap=cmap(numEns,:);
% cmap(end+1,:)=rgb('grey');
if numEns==1
    cmap{1} = rgb('Grey');
else
cmap = colorMapPicker(numEns,outVars.defaultColorMap);
end
cmap{end+1} = rgb('Firebrick');

clear p

p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
% bar(x, avg)
% hold on
% er = errorbar(x, avg, err);
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% hold off
ylabel('Population Response')
% xticklabels(uniqueEns)
% xticks = 1:6;
title('Mean Population Response to Ensemble')
% xlabel('Ensemble Size')
set(gcf(),'Name','Mean population response to holo')
% ns

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 2;
p{2}(2).LineWidth = 2;
% p{2}(2).

r = refline(0);
r.LineStyle=':';
r.Color = rgb('grey');
r.LineWidth = 1;

pValEnselbeSize = anovan(popResponseEns(ensemblesToUse),numCellsEachEns(ensemblesToUse)','display','off');

disp(['Anova pVal between sizes: ' num2str(pValEnselbeSize)]);
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))
uniqueEns(end+1) = 0; 
for i=1:size(data,2)
    %     prs = ranksum(data{i},0);
    psr = signrank(data{i});
    
    [h pp ] = ttest(data{i},0);
    disp(['Size: ' num2str(uniqueEns(i)) ', n=' num2str(numel(data{i}))...
        '. Mean: ' num2str(nanmean(data{i}),2) ' ' char(177) ' ' num2str(sem(data{i}),2)...
        '. Signed Rank: ' num2str(psr,2) '. ttest: ' num2str(pp,2)])
end