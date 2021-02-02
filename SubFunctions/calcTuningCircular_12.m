function [All, outVars] = calcTuningCircular_12(All, outVars)

% settings
numExps = numel(All);

if ~isfield(All(1).out.anal, 'oriCurve')
    error("No ori curve, can't calculate OSI.")
end

for ind = 1:numExps
    [tune_result, var_result, curves_result] = tuningByCircular_12(All(ind).out.anal.oriCurve);
    All(ind).out.anal.circTuning = tune_result;
    All(ind).out.anal.circVar = var_result;
    All(ind).out.anal.circCurves = curves_result;
    tunings2save{ind} = tune_result;
    var2save{ind} = var_result;
    curve2save{ind} = curves_result;
end

outVars.circTuning = tunings2save;
outVars.circVar = var2save;
outVars.circCurves = curve2save;