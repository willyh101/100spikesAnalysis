numExpts = numel(All);
c=0;
clear mRespNoVis mRespVis
for ind = 1:numExpts
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];

    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    cellList = 1:numel(ROIinArtifact);
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    pVisR = All(ind).out.anal.pVisR;
    trialsToUse = All(ind).out.anal.defaultTrialsToUse;
    
    
    
    for i = 1:numel(us)
        u = us(i);
        
        h = All(ind).out.exp.stimParams.roi{i};
        if h>0
            tg = All(ind).out.exp.holoTargets{h};
            tg(isnan(tg))=[];
        else
            tg=[];
        end
        
        c=c+1;
        
        if i==1
            cellsToUse = ~ROIinArtifact' & pVisR < 0.05;
        else
            cellsToUse = ~ROIinArtifact' & pVisR < 0.05 & ~offTargetRisk(h,:) & ~ismember(cellList,tg);
        end
            
        
        % no vis case
        v = 1;
        respsubt1 = squeeze(respMat(i,v,cellsToUse)) - squeeze(baseMat(i,v,cellsToUse));
        mRespNoVis(c) = squeeze(mean(respsubt1));
        
        % any vis case
        v = vs>1;
        respsubt2 = squeeze(mean(respMat(i,v,cellsToUse))) - squeeze(mean(baseMat(i,v,cellsToUse)));
        mRespVis(c) = squeeze(mean(respsubt2));
    end
end
        
        
        