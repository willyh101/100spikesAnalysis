function [All, outVars] = compareRedCellsTuning(All, outVars)

clear coTunedRedResp notCoTunedRedResp coTunedRedEnsResp notCoTunedRedEnsResp

mRespRed = outVars.mRespRed; % {expt}(cells, ensemble)
umouse = outVars.uMouse; % which mouse/expt each ensemble belongs to
ensTunings = rmmissing(unique(outVars.ensPO)); % must have > Matlab2018-ish to nix nans this way

% maybe make a table so it's easier?

% find ensembles per mouse
mice = diff([0, find(diff(umouse)), numel(umouse)]);

% make mRespEns a cell array by animal so it's easier to index
mRespRed = cellfun(@cell2mat, mat2cell(mRespRed, 1, mice), 'un', 0);

% for each tuning what is the mean resp for red cell that are tuned and not
% tuned to that ensStim
for i = 1:numel(ensTunings)
    tuning = ensTunings(i); % is in degrees

    for ind = 1:numel(All)
        
        % red cell tuning is an arr of tunings, the length of redVisCells,
        % but we want their positional index out of all red cells so we can
        % then index into mRespRed
        redCells = find(All(ind).out.red.isRed); % all red cells logical arr (used to get index)
        redVisCells = All(ind).out.red.isVisCells; % red and vis responsive cells, logical arr
        redVisIdx = find(ismember(redCells, redVisCells));
        redCellTuning = All(ind).out.red.redTuningOri;
        
        % nans will show up as not-co-tuned but are visually responsive
        % but they won't appear as a 'tuning' because they are removed as a
        % uniquely shot tuning
        isCoTuned = redVisIdx(ismember(redCellTuning, tuning));
        isNotCoTuned = redVisIdx(~ismember(redCellTuning, tuning));
        
        % so now, get the co-tuned, and not-co-tuned mean responses
        coTunedRedResp{ind} = mean(mRespRed{ind}(isCoTuned, :), 1);
        notCoTunedRedResp{ind} = mean(mRespRed{ind}(isNotCoTuned, :), 1);
      
    end
    
%     % mean across mice to get mean response
%     coTunedRedEnsResp(:, i) = mean(coTunedRedResp);
%     notCoTunedRedEnsResp(:, i) = mean(notCoTunedRedResp);
        coTunedRedEnsResp(:, i) = cellfun(@mean, coTunedRedResp);
        notCoTunedRedEnsResp(:, i) = arrayfun(@mean, notCoTunedRedResp{:});
    
end
    
outVars.coTunedRedEnsResp = coTunedRedEnsResp;
outVars.notCoTunedRedEnsResp = notCoTunedRedEnsResp;


