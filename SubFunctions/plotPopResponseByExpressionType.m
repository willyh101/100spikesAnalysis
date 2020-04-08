function plotPopResponseByExpressionType(All,outVars);
%%
sepSizeColors=0;

ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
popResponseEns = outVars.popResponseEns;
noStimPopResp = outVars.noStimPopResp;
ensIndNumber = outVars.ensIndNumber;
ensExpressionType = outVars.ensExpressionType;

uniqueExpressionTypes = outVars.uniqueExpressionTypes;

f = figure(8);
clf

ExpTypes = unique(ensExpressionType(ensemblesToUse));

IndsUsed = outVars.IndsUsed;

names=[];
for Ind = 1:numel(All)
    names{Ind}=strrep(All(Ind).out.info.mouse, '_', '.');
end


[ensSizes, ~, numCellsIdx] = unique(numCellsEachEns(ensemblesToUse));
colorList = colorMapPicker(numel(ensSizes),'parula');

for i = 1:numel(ExpTypes);
    E = ExpTypes(i);
    
    sp(i) = subplot(1,numel(ExpTypes),i);
    
    datToUse = popResponseEns(ensemblesToUse & ensExpressionType==E);
    ensIDsToUse = ensIndNumber(ensemblesToUse & ensExpressionType==E);
    [uniqueToUse, ~, ensIDsToUseTemp] = unique(ensIDsToUse);
    numCellsEachEns(ensemblesToUse & ensExpressionType==E);
    numCellsToUse = numCellsIdx(ensExpressionType(ensemblesToUse)==E);
    %     [~, ~, numCellsToUseTemp] = unique(numCellsToUse);
    
    if sepSizeColors
    cmapToUse = cell2mat(colorList(numCellsToUse)');
    scatter(ensIDsToUseTemp,datToUse,[],cmapToUse,'filled')
    else
    scatter(ensIDsToUseTemp,datToUse,[],rgb('grey'),'filled')
    end
    
    hold on
    r= refline(0);
    r.LineStyle = ':';
    r.LineWidth = 2;
    r.Color = rgb('Grey');
    
    popRespMean=[];
    popRespSEM=[];
    for k=1:numel(uniqueToUse)
        popRespMean(k) = nanmean(datToUse(ensIDsToUseTemp==k));
        popRespSEM(k) = sem(datToUse(ensIDsToUseTemp==k))*1.96;
    end
    
    
    e = errorbar(1:numel(uniqueToUse),popRespMean,popRespSEM);
    e.LineStyle='none';
    e.LineWidth =2;
    e.Color = rgb('red');
    
    xticks(1:numel(uniqueToUse))
    xticklabels(names(uniqueToUse));
    xtickangle(45)
    
    title(uniqueExpressionTypes{E})
end
linkaxes(sp,'y')