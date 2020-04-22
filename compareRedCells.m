function [All, outVars] = compareRedCells(All, outVars, opts)

% opts.redCellCompare = 'isVis';
cond = opts.redCellCompare;

ensemblesToUse = outVars.ensemblesToUse;
mRespRed = outVars.mRespRed;
mRespNotRed = outVars.mRespNotRed;
umouse = outVars.uMouse;

mice = diff([0, find(diff(umouse)), numel(umouse)]);

mRespRed = cellfun(@cell2mat, mat2cell(mRespRed, 1, mice), 'un', 0);
mRespNotRed = cellfun(@cell2mat, mat2cell(mRespNotRed, 1, mice), 'un', 0);

meanRedCond = [];
semRedCond = [];
meanRedInv = [];
semRedInv = [];

for ind = 1:numel(All)
    redCells = find(All(ind).out.red.isRed);
    redCondCells = All(ind).out.red.(cond);
    
    assert(isa(redCondCells, 'logical'), 'Input array must be a list of cells. Do not use logical indexing.')
    
    redCondLogical = zeros(1, numel(redCells));
    redCondLogical(redCondCells) = 1;
    
    redNotVisCells = find(redCondLogical & All(ind).out.red.isRed);
    [~, redCond] = ismember(redCondCells, redCells); 
    
    meanRedVis = [meanRedVis mean(mRespRed{ind}(redCond, :))];
    semRedVis = [semRedVis sem2(mRespRed{ind}(redCond, :))];
    
    meanRedInv = [meanRedInv mean(mRespRed{ind}(redCond, :))];
    semRedInv = [semRedInv sem2(mRespRed{ind}(redCond, :))];
    
end


%%
ensemblesToUse = outVars.ensemblesToUse;
mRespRed = outVars.mRespRed;
mRespNotRed = outVars.mRespNotRed;

cond = 'isVis'; % must be logical

for ind = numel(All)
    cond = All(ind).out.red.(cond);
    
    redCells = find(All(ind).out.red.isRed);
%     redVisCells = All(ind).out.red.isVisCells;
    % more generically, just...
    condCells = All(ind).out.red.(cond);
    
    assert(isa(condCells, 'logical'), 'Input array must be a list of cells. Do not use logical indexing.')
    
    invCondRedCells = find(~condCells & All(ind).out.red.isRed);
    
%     redNotVisCells = find(~All(ind).out.red.isVis & All(ind).out.red.isRed);
    [~, redCond] = ismember(redVisCells, redCells);
    
    
    
end
