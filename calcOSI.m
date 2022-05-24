function [All, outVars] = calcOSI(All, outVars)

% settings
numExps = numel(All);

assert(isfield(All(1).out.anal, 'oriCurve'), ...
    "No ori curve, can't calculate OSI.")

for ind = 1:numExps
    curve = All(ind).out.anal.oriCurve;
    [~, po] = max(curve, [], 1);
    curve(:, po==1) = nan;
    curve(1,:) = [];
    result = osi(curve);
    All(ind).out.anal.osi = result.osi;
    All(ind).out.anal.poFromOSI = result.po;
    osi2save{ind} = result.osi;
    po2save{ind} = result.po;
end

outVars.osi = osi2save;
outVars.poFromOSI = po2save;