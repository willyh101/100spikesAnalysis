function [All, outVars] = calcOSI(All, outVars)

% settings
numExps = numel(All);

if ~isfield(All(1).out.anal, 'oriCurve')
    error("No ori curve, can't calculate OSI.")
end

for ind = 1:numExps
    result = osi(All(ind).out.anal.oriCurve);
    All(ind).out.anal.osi = result;
    osi2save{ind} = result;
end

outVars.osi = osi2save;
