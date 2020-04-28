function [All, outVars] = compareAllCellsTuning(All, outVars, opts)

visAlpha = opts.visAlpha;

mRespEns = outVars.mRespEns; % ensemble mResp group by experiment/outfile with indiv cells
ensExpID = outVars.uExptsEns;
ensTunings = outVars.ensPO; % PO of all ensembles
uEnsTunings = rmmissing(unique(ensTunings)); % get rid of nans since they are all unique

% transforms/groups ensemble vectors into experiment-group cell array
ensExpID = diff([0, find(diff(ensExpID)), numel(ensExpID)]);
ensTunings = arrayfun(@cell2mat, mat2cell(ensTunings, 1, ensExpID), 'un', 0);
mRespEns = cellfun(@cell2mat, mat2cell(mRespEns, 1, ensExpID), 'un', 0);


for ind = 1:numel(All)
    isVis = All(ind).out.anal.pVisR < visAlpha;
    try
        isTuned = All(ind).out.anal.pVisT < visAlpha;
    catch
        isTuned = isVis;
    end
    isVisCells = find(isVis);
    isTunedCells = find(isTuned);
    pref = idx2ori(All(ind).out.anal.prefOri(isVisCells), [nan 0:45:315]);
    orthos = idx2ori(All(ind).out.anal.orthoOri(:, isVisCells), [nan 0:45:315]);
    
    clear allResp isoResp orthoResp notTunedREsp notVisResp
    for tune = 1:numel(uEnsTunings)
        tuning = uEnsTunings(tune);
        ensTuned = ensTunings{ind} == tuning;
        
        isIso = isVisCells(ismember(pref, tuning));
        isOrtho = isVisCells(any(ismember(orthos, tuning)));
        notCoTuned = isVisCells(~ismember(pref, tuning)); % should be excluding cells in the not tuned condition, I think?
%         notTuned = ~isTuned;
        notVis = ~isVis;
        
        allResp{ind}(tune, :) = nanmean(mRespEns{ind}, 1);
        allResp{ind}(tune, ~ensTuned) = nan;
        
        isoResp{ind}(tune, :) = nanmean(mRespEns{ind}(isIso, :), 1);
        isoResp{ind}(tune, ~ensTuned) = nan;
        
        orthoResp{ind}(tune, :) = nanmean(mRespEns{ind}(isOrtho, :), 1);
        orthoResp{ind}(tune, ~ensTuned) = nan;
        
        notCoTunedResp{ind}(tune, :) = nanmean(mRespEns{ind}(notCoTuned, :), 1);
        notCoTunedResp{ind}(tune, ~ensTuned) = nan;
        
%         notTunedResp{ind}(tune, :) = nanmean(mRespEns{ind}(notTuned, :), 1);
%         notTunedResp{ind}(tune, ~ensTuned) = nan;
        
        notVisResp{ind}(tune, :) = nanmean(mRespEns{ind}(notVis, :), 1);
        notVisResp{ind}(tune, ~ensTuned) = nan;
        
    end
    
    allEnsResp{ind} = nanmean(allResp{ind}, 1);
    isoEnsResp{ind} = nanmean(isoResp{ind}, 1);
    orthoEnsResp{ind} = nanmean(orthoResp{ind}, 1);
    notCoTunedEnsResp{ind} = nanmean(notCoTunedResp{ind}, 1);
%     notTunedEnsResp{ind} = nanmean(notTunedResp{ind}, 1);
    notVisEnsResp{ind} = nanmean(notVisResp{ind}, 1);
    
end

outVars.allEnsResp = cell2mat(allEnsResp);
outVars.isoEnsResp = cell2mat(isoEnsResp);
outVars.orthoEnsResp = cell2mat(orthoEnsResp);
outVars.notCoTunedResp = cell2mat(notCoTunedEnsResp);
% outVars.notTunedEnsResp = cell2mat(notTunedEnsResp);
outVars.notVisEnsResp = cell2mat(notVisEnsResp);
