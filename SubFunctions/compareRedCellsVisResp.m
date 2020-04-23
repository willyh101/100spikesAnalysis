function [All, outVars] = compareRedCellsVisResp(All, outVars)

ensemblesToUse = outVars.ensemblesToUse;

mRespRed = outVars.mRespRed;
mRespNotRed = outVars.mRespNotRed;
umouse = outVars.uMouse;

mice = diff([0, find(diff(umouse)), numel(umouse)]);

mRespRed = cellfun(@cell2mat, mat2cell(mRespRed, 1, mice), 'un', 0);
mRespNotRed = cellfun(@cell2mat, mat2cell(mRespNotRed, 1, mice), 'un', 0);

meanRedVis = [];
semRedVis = [];
meanRedNotVis = [];
semRedNotVis = [];

for ind = 1:numel(All)
    
    redCells = find(All(ind).out.red.isRed);
    redVisCells = All(ind).out.red.isVisCells; 
    
    redVisIdx = ismember(redCells, redVisCells);
    redNotVisIdx = ~ismember(redCells, redVisCells);
    
    meanRedVis = [meanRedVis mean(mRespRed{ind}(redVisIdx, :), 1)];    
    meanRedNotVis = [meanRedNotVis mean(mRespRed{ind}(redNotVisIdx, :), 1)];
    
end

outVars.popRespEnsRedVisResp = meanRedVis;
outVars.popRespEnsRedNotVisResp = meanRedNotVis;


