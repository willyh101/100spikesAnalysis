function result = osi(curve)
% do not include catch condition!
% should be condition x cell
% 
% try
%     assert( size(curve,1) < size(curve,2))
% catch
%     warning('Dim 1 of tuning curve is not greater than dim 2 of tuning curve. Make sure tuning curve is ori x cell.')
% end
    
curve = curve - min(curve);
curve = makeTuningCurve180(curve);

[~, prefOri]= max(curve, [], 'omit');

orthoOri = mod(prefOri + 2, 4);
orthoOri(orthoOri==0) = 4;

OSI = zeros(size(prefOri));

for i=1:numel(prefOri)
    OSI(i) = (curve(prefOri(i),i) - curve(orthoOri(i),i)) ...
        / (curve(prefOri(i),i) + curve(orthoOri(i),i));
end

result = OSI;