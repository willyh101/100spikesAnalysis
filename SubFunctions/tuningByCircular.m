function [PO, CV] = tuningByCircular(oris, data)

% remove the cell if its max is the catch
% I think that keeping it in messes things up... (ie. not a circle)
data = data'; % needs to be cells x oris
[~,p] = max(data,[],2);
data(p==1,:)=[]; % remove cellsif match is grey screen
data(:,1)=[]; % remove the data from array
data = data - min(data); % remove negative values

% convert to rads
ori = circ_axial(circ_ang2rad(oris), 2);

% calculate bin spacing in rads
dori = diff(ori(1:2));

% calculations
ncells = size(data,1);
for c = 1:ncells
    d = data(c,:);
    % get mean and circular variance
    m(c) = circ_mean(ori, d, 2);
    v(c) = circ_var(ori, d, dori, 2);
    s(c) = circ_std(ori, d, dori, 2);
    rt(c) = circ_rtest(ori,d, dori); %idk if this is working
%     ot(c) = circ_otest(ori, [], d); % needs integers
end

PO = circ_rad2ang(circ_axial(m));
CV = circ_rad2ang(circ_axial(v));