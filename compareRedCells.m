function [All, outVars] = compareRedCells(All, outVars, opts)

% opts.redCellCompare = 'isVis';
cond = opts.redCellCompare;

mRespRed = outVars.mRespRed;
mRespNotRed = outVars.mRespNotRed;
umouse = outVars.uMouse;

mice = diff([0, find(diff(umouse)), numel(umouse)]);

mRespRed = cellfun(@cell2mat, mat2cell(mRespRed, 1, mice), 'un', 0);
mRespNotRed = cellfun(@cell2mat, mat2cell(mRespNotRed, 1, mice), 'un', 0);

meanRedCond = [];
meanRedInv = [];

for ind = 1:numel(All)
    redCells = find(All(ind).out.red.isRed);
    
    if isa(All(ind).out.red.(cond), 'logical')
        redCondCells = find(All(ind).out.red.(cond));
    elseif isa(All(ind).out.red.(cond), 'double')
        redCondCells = All(ind).out.red.(cond);
    else
        error('Input array must be of type logical or double (list of cells of logical index of cells).')
    end
        
    redCond = ismember(redCells, redCondCells); 
    redInv = ~ismember(redCells, redCondCells);
    
    meanRedCond = [meanRedVis mean(mRespRed{ind}(redCond, :))];
    meanRedInv = [meanRedInv mean(mRespRed{ind}(redCond, :))];
    
end

result = ['popRespEnsRedBy_' cond];
resultInv = ['popRespEnsRedBy_Inv_' cond];

outVars.(result) = meanRedCond;
outVars.(resultInv) = meanRedInv;
