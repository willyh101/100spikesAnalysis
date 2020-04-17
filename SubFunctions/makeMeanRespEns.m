function outVars = makeMeanRespEns(All, outVars)

% get mResp of red cells

numExps = numel(All);

c = 0;
for ind = 1:numExps
    
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    
    us = unique(All(ind).out.exp.stimID);
    
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;
    
    isRed = All(ind).out.red.isRed;
    
    v = 1;
    
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    
    for i = 1:numStims
        
        holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
        c=c+1;
        
        mRespRed{c} = squeeze(respMat(i, v, :) - ...
           baseMat(i, v, :));
        mRespNotRed{c} = squeeze(respMat(i, v, :) - ...
            baseMat(i, v, :));
        
        if i==1
            mRespRed{c}(~ROIinArtifact) = nan;
            mRespNotRed{c}(~ROIinArtifact) = nan;
        else
            mRespRed{c}(~offTargetRisk(holo, :) & ~ROIinArtifact) = nan;
            mRespNotRed{c}(~offTargetRisk(holo, :) & ~ROIinArtifact) = nan;
        end
        
        mRespRed{c} = mRespRed{c}(isRed);
        mRespNotRed{c} = mRespNotRed{c}(~isRed);
        
        isControl(c) = i==1;
        
        umouse(c) = ind;
    end
    
end

% a = mRespRed


outVars.mRespRed = mRespRed(~isControl);
outVars.mRespNotRed = mRespNotRed(~isControl);
outVars.uMouse =  umouse;