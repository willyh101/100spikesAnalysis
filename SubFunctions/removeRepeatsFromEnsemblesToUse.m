function [ensemblesToUse, outVars] = removeRepeatsFromEnsemblesToUse(ensemblesToUse,outVars)
subRepeats = outVars.repeatIdentifier(ensemblesToUse);
ensToUseList = find(ensemblesToUse);
[a,b]=histc(subRepeats,unique(subRepeats));
removedRepeats = zeros(size(ensemblesToUse));
for i = 2:numel(a);
    if a(i)>1
        problemCells = ensToUseList(b==i);
        
        if sum(outVars.hzEachEns(problemCells)==10)==1;
        removedRepeats(problemCells(outVars.hzEachEns(problemCells)~=10))=1;
        else
            removedRepeats(problemCells(2:end))=1;
        end
    end
end
ensemblesToUse = ensemblesToUse & ~removedRepeats;
outVars.removedRepeats = removedRepeats;