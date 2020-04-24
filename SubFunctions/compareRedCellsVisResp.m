function [All, outVars] = compareRedCellsVisResp(All, outVars, opts)

mRespRed = outVars.mRespRed;
mRespNotRed = outVars.mRespNotRed;
umouse = outVars.uMouse;

mice = diff([0, find(diff(umouse)), numel(umouse)]);

mRespRed = cellfun(@cell2mat, mat2cell(mRespRed, 1, mice), 'un', 0);
mRespNotRed = cellfun(@cell2mat, mat2cell(mRespNotRed, 1, mice), 'un', 0);

meanRedVis = [];
meanRedNotVis = [];
meanOtherVis = [];
meanOtherNotVis = [];

for ind = 1:numel(All)
    
    % red cells
    redCells = find(All(ind).out.red.isRed);
    redVisCells = All(ind).out.red.isVisCells; 
    
    redVisIdx = ismember(redCells, redVisCells);
    redNotVisIdx = ~ismember(redCells, redVisCells);
    
    meanRedVis = [meanRedVis mean(mRespRed{ind}(redVisIdx, :), 1)];    
    meanRedNotVis = [meanRedNotVis mean(mRespRed{ind}(redNotVisIdx, :), 1)];
    
    % other cells
    otherCells = find(~All(ind).out.red.isRed);
    otherVisCells = All(ind).out.red.isVisOtherCells;
    
    otherVisIdx = ismember(otherCells, otherVisCells);
    otherNotVisIdx = ~ismember(otherCells, otherVisCells);
    
    meanOtherVis = [meanOtherVis, mean(mRespNotRed{ind}(otherVisIdx, :), 1)];
    meanOtherNotVis = [meanOtherNotVis, mean(mRespNotRed{ind}(otherNotVisIdx, :), 1)];
    
end

outVars.popRespEnsRedVisResp = meanRedVis;
outVars.popRespEnsRedNotVisResp = meanRedNotVis;
outVars.popRespEnsOtherVisResp = meanOtherVis;
outVars.popRespEnsOtherNotVisResp = meanOtherNotVis;


