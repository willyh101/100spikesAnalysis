us = unique(All(ind).out.exp.stimID);
for k=1:numTrials
    s = All(ind).out.exp.stimID(k);
    sidx = find(us==s);
    
    % h is the hologram ID
    h = All(ind).out.exp.stimParams.roi{sidx};
    % htg are the matched cells shot during that stimID;
    % rtg are the rois (index of stimCoM / stimDepth);
    if h==0
        htg=[];
    else
        htg = All(ind).out.exp.holoTargets{h};
        rtg = All(ind).out.exp.rois{h};
        htg(isnan(htg))=[];
    end
end


%% the right way to ID holos / holoTargets

for i  = 1:ensToUse
    e = ensToUse(i);
    ind = ensIndNumber(e);
    
    s = outVars.ensHNumber(e); %i think this is the order of resp
    h = All(ind).out.exp.stimParams.roi{s};
    htg = All(ind).out.exp.holoTargets{h};
    htg(isnan(htg))=[];
    
    rtg = All(ind).out.exp.rois{h};
    
    us = unique(All(ind).out.exp.stimID);

    sID = us(s); %StimID value for this ensemble
end