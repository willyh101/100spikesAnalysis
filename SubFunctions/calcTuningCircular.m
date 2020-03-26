function [All, outVars] = calcTuningCircular(All, outVars)

% settings
numExps = numel(All);

if ~isfield(All(1).out.anal, 'oriCurve')
    error("No ori curve, can't calculate OSI.")
end

oris = 0:45:315;
for ind = 1:numExps
    [tune_result, var_result] = tuningByCircular(oris, All(ind).out.anal.oriCurve);
    All(ind).out.anal.circTuning = tune_result;
    All(ind).out.anal.circVar = var_result;
    tunings2save{ind} = tune_result;
    var2save{ind} = var_result;
end

outVars.circTuning = tunings2save;
outVars.circVar = var2save;