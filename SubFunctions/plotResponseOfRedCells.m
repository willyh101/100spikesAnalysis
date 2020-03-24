function plotResponseOfRedCells(All,outVars,opts)

numExps = numel(All);


for ind = 1:numExps;
    
    
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);
    
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    
    isRed = All(ind).out.red.isRed;
    
    v=1;
        
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    clear popRespRed popRespNotRed popRespRedSEM popRespNotRedSEM
 for i= 1:numStims
            holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
            if i==1;
                cellsToUse = ~ROIinArtifact';
            else
                cellsToUse = ~ROIinArtifact'  & ~offTargetRisk(holo,:);
            end
            popRespRed(i) = nanmean(squeeze(respMat(i,v,cellsToUse & isRed) - baseMat(i,v,cellsToUse & isRed)));
            popRespRedSEM(i) = sem(squeeze(respMat(i,v,cellsToUse & isRed) - baseMat(i,v,cellsToUse & isRed)));
            popRespNotRed(i) = mean(squeeze(respMat(i,v,cellsToUse & ~isRed) - baseMat(i,v,cellsToUse & ~isRed)));
            popRespNotRedSEM(i) = sem(squeeze(respMat(i,v,cellsToUse & ~isRed) - baseMat(i,v,cellsToUse & ~isRed)));

 end
 
%  figure(23);clf
%  plot(popRespRed,'o')
%  hold on
%  plot(popRespNotRed,'o')
% 
%  plot(popRespRed,popRespNotRed,'o')
 
end