function [outVars] = identifyDuplicateHolos(All,outVars)
disp('Identifying Repeated Holograms')
numExps = numel(All);

c=0;

repeatIdentifier = [];
for ind = 1:numExps
    us = unique(All(ind).out.exp.stimID);
    thisExpROIs=All(ind).out.exp.rois;
    
    for k=1:numel(us)
        s = us(k);
        sidx = find(us==s);
        
        % h is the hologram ID
        h = All(ind).out.exp.stimParams.roi{sidx};
        % htg are the matched cells shot during that stimID;
        % rtg are the rois (index of stimCoM / stimDepth);
        if h==0
            rtg=[];
        else
            c=c+1;
            rtg = All(ind).out.exp.rois{h};
            repeatsList =cellfun(@(x) isequal(unique(x(:)),unique(rtg(:))),thisExpROIs,'uniformOutput',1);
            numReps = sum(repeatsList);
            if numReps==1;
                repeatIdentifier(c)=0;
            else
                repeatIdentifier(c)=ind*1000+find(repeatsList,1); %uniqueIdentifier 
            end
        end
    end
end
outVars.repeatIdentifier = repeatIdentifier; 