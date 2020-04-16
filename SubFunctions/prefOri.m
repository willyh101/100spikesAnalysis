function [preferred, orthogonal] = prefOri(curve)
% calculates the preferred and orthogonal orientations
% tuning curve can include catch condition
% returns indices to get degrees for each

assert( size(curve,1) < size(curve,2), ... 
    ['Dim 1 of tuning curve must be greater than dim 2 of tuning curve. '...
    'Make sure tuning curve is ori x cell.'] )

conds = size(curve, 1);
ortho_offset = floor(conds / 2);

[~, pref] = max(curve);

% mark the catch conditions with nan
pref(pref == 1) = nan;

ortho1 = mod(pref + ortho_offset, conds);
ortho2 = mod(pref - ortho_offset, conds);

% fix for making a 0 index...
ortho1(ortho1 == 0) = conds;
ortho2(ortho2 == 0) = conds;

% make nans 1s again
pref(isnan(pref)) = 1;
ortho1(isnan(ortho1)) = 1;
ortho2(isnan(ortho2)) = 1;

% return
preferred = pref;
orthogonal = [ortho1; ortho2];