function [outVars] = detectShotRedCells(All,outVars)
numExps = numel(All);

c=0;
for ind = 1:numExps
    
    us = unique(All(ind).out.exp.stimID);
    
    for i = 1:numel(us)
%         s = us(i);
        
        h = All(ind).out.exp.stimParams.roi{i};
        if h==0
            htg=[];
            
        else
            htg = All(ind).out.exp.holoTargets{h};
            htg(isnan(htg))=[];
            c=c+1;
            if isfield(All(ind).out, 'red')
                ensHasRed(c) =any(All(ind).out.red.isRed(htg));
            else
                ensHasRed(c)=0;
            end
            
        end
    end
end
outVars.ensHasRed = ensHasRed; 