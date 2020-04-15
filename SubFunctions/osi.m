function result = osi(curve)
% do not include catch condition!
% should be condition x cell

assert( size(curve,1) < size(curve,2), ... 
    ['Dim 1 of tuning curve must be greater than dim 2 of tuning curve. '...
    'Make sure tuning curve is ori x cell.'] )
    
curve = curve - min(curve);
curve = makeTuningCurve180(curve);

[~, prefOri]= max(curve);

orthoOri = mod(prefOri + 2, 4);
orthoOri(orthoOri==0) = 4;

OSI = zeros(size(prefOri));

for i=1:numel(prefOri)
    OSI(i) = (curve(prefOri(i),i) - curve(orthoOri(i),i)) ...
        / (curve(prefOri(i),i) + curve(orthoOri(i),i));
end

result = OSI;