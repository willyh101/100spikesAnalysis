function [outVars] = missedTargetDetector(All,outVars,opts)
numExps = numel(All);

c=0;
ensMissedTargetF = [];
ensMissedTarget =[];
for ind = 1:numExps
    us = unique(All(ind).out.exp.stimID);
    for k=1:numel(us)
        s = us(k);
        sidx = find(us==s);
        
        % h is the hologram ID
        h = All(ind).out.exp.stimParams.roi{sidx};
        % htg are the matched cells shot during that stimID;
        % rtg are the rois (index of stimCoM / stimDepth);
        if h==0
            htg=[];
            rtg=[];
        else
            c=c+1;
            htg = All(ind).out.exp.holoTargets{h};
            rtg = All(ind).out.exp.rois{h};
            htg(isnan(htg))=[];
            
            F = 1- (numel(htg)/numel(rtg));
            ensMissedTargetF(c) = F; %Fraction of targets not detected by s2p
            ensMissedTarget(c)  = F > opts.FractionMissable;
            
        end
    end
    
    
end

outVars.ensMissedTargetF    = ensMissedTargetF;
outVars.ensMissedTarget     = ensMissedTarget;

disp([num2str(sum(ensMissedTarget)) ' Ensembles Failed. ' num2str(mean(ensMissedTarget)*100) ' %'])