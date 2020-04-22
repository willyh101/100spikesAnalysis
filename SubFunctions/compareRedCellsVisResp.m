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
    
%     redNotVisCells = find(~All(ind).out.red.isVis & All(ind).out.red.isRed);
    redVisIdx = ismember(redCells, redVisCells);
%     [~, redVisIdx] = ismember(redVisCells, redCells);
    redNotVisIdx = ~ismember(redCells, redVisCells);
    
    meanRedVis = [meanRedVis mean(mRespRed{ind}(redVisIdx, :), 1)];
%     semRedVis = [semRedVis sem2(mRespRed{ind}(redVisIdx, :))];
    
    meanRedNotVis = [meanRedNotVis mean(mRespRed{ind}(redNotVisIdx, :), 1)];
%     semRedNotVis = [semRedNotVis sem2(mRespRed{ind}(redNotVisIdx, :))];
end

outVars.popRespEnsRedVisResp = meanRedVis;
outVars.popRespEnsRedNotVisResp = meanRedNotVis;


