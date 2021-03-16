function [All, outVars] = calcOSI_GH(All, outVars)

% settings
numExps = numel(All);

assert(isfield(All(1).out.anal, 'oriCurve'), ...
    "No ori curve, can't calculate OSI.")

for ind = 1:numExps
    curve = All(ind).out.anal.oriCurve;
    [~, po] = max(curve, [], 1);
    curve(:, po==1) = nan;
    curve(1,:) = [];
    [result, GOSI] = osi_GH(curve);
    All(ind).out.anal.osi = result;
    All(ind).out.anal.gosi = GOSI;
    osi2save{ind} = result;
    gosi2save{ind} = GOSI;
end

outVars.osi = osi2save;
outVars.gosi = gosi2save;
