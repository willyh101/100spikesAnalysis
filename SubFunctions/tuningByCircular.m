% function [PO, CV] = tuningByCircular(oris, data)

oris = 0:45:315;
data = outVars.tuningCurves{1}';
% data = data'; % needs to be cells x oris
data(:,1) = []; % remove the catch condition
% data = data - min(data);

% convert to rads
ori = circ_axial(circ_ang2rad(oris), 2);
% ori = circ_axial(oris, 2);

% calculate bin spacing in rads
dori = diff(ori(1:2));

% calculations
ncells = size(data,1);
for c = 1:ncells
    d = data(c,:);
    
    m(c) = circ_mean(ori, d, 2);
    v(c) = circ_var(ori, d, dori, 2);
    
    mw = max(d);
    r(c,:) = circ_r(ori, d, dori) * mw;
    phi(c,:) = circ_mean(ori,d);
    
    % tests might not work for non-ordinal data
    rt(c) = circ_rtest(ori,d, dori); %idk if this is working
%     ot(c) = circ_otest(ori, [], d); % needs integers
end

PO = circ_rad2ang(circ_axial(m));
% PO = circ_rad2ang(m);
CV = circ_rad2ang(circ_axial(v));
% CV = circ_rad2ang(v);
% CV = v;

figure(1)
clf
subplot(1,2,1)
histogram(PO,100)
title('Circular Mean')
xlabel('Orientation')
ylabel('Count')
subplot(1,2,2)
histogram(CV,100)
title('Circular Variance')
xlabel('Degrees')