% opts.redCellCompare = 'isVis';
cond = opts.redCellCompare;

ensemblesToUse = outVars.ensemblesToUse;
mRespRed = outVars.mRespRed;
mRespNotRed = outVars.mRespNotRed;
umouse = outVars.uMouse;

mice = diff([0, find(diff(umouse)), numel(umouse)]);

mRespRed = cellfun(@cell2mat, mat2cell(mRespRed, 1, mice), 'un', 0);
mRespNotRed = cellfun(@cell2mat, mat2cell(mRespNotRed, 1, mice), 'un', 0);

meanRedVis = [];
semRedVis = [];

for ind = 1:numel(All)
    redCells = find(All(ind).out.red.isRed);
    redCondCells = All(ind).out.red.(cond);
    
    assert(~isa(redCondCells, 'logical'), 'Input array must be a list of cells. Do not use logical indexing.')
    
    redCondLogical = zeros(1, numel(redCells))
    redCondLogical(redCondCells) = 1;
    
    redNotVisCells = find(~All(ind).out.red.isVis & All(ind).out.red.isRed);
    [~, redVisIdx] = ismember(redVisCells, redCells);
    
    meanRedVis = [meanRedVis mean(mRespRed{ind}(redVisIdx, :))];
    semRedVis = [semRedVis sem2(mRespRed{ind}(redVisIdx, :))];
end
