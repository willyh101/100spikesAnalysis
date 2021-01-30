function [All, outVars] = calcOSI_12(All, outVars)

% settings
numExps = numel(All);

assert(isfield(All(1).out.anal, 'oriCurve'), ...
    "No ori curve, can't calculate OSI.")

for ind = 1:numExps
    curve = All(ind).out.anal.oriCurve;
    [~, po] = max(curve, [], 1);
    curve(:, po==1) = nan;
    curve(1,:) = [];
    result = osi_12(curve);
    All(ind).out.anal.osi = result;
    osi2save{ind} = result;
end

outVars.osi = osi2save;
