function plotVisPercents(All,outVars)
%%
numExps = numel(All);

visThreshold = 0.05;
indExpressionType = outVars.indExpressionType; 
[a indOrderToUse] = sort(indExpressionType);

names=[];
for Ind = 1:numel(All)
    names{Ind}=strrep(All(Ind).out.info.mouse, '_', '.');
end

c=0;
clear overalVis targetVis shotVis
for ind = indOrderToUse
   c=c+1;
    pVisR = All(ind).out.anal.pVisR;
    
    targetedCells = All(ind).out.exp.targetedCells;
    targetedCells(isnan(targetedCells))=[];
    
    shotCells = unique([All(ind).out.exp.holoTargets{:}]);
    shotCells(isnan(shotCells))=[];
    
    
    overalVis(c) = mean(pVisR<visThreshold);
    targetVis(c) = mean(pVisR(targetedCells)<visThreshold);
    shotVis(c) = mean(pVisR(shotCells)<visThreshold);
    
end

%% Plot
toPlot =[overalVis; targetVis; shotVis]';

figure(43);clf
bar(toPlot)

xticks(1:c)

reorderUnique = outVars.uniqueExpressionTypes(a);
reordernames = names(indOrderToUse)';

for i = 1:numel(reorderUnique)
    labelsToPrint{i} = [reorderUnique{i} '\newline     ' reordernames{i}];
end
xticklabels(labelsToPrint)
xtickangle(90)
legend('All Cells','TargetedCells','Shot cells')

ylabel('Proportion of Cells Vis Responsive')
xlabel('Expression Strategy')