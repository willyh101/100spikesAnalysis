function outVars = makeMeanRespEnsByCell(All, outVars)

numExps = numel(All);

c=0;c2=0;
for ind = 1:numExps
    
    us = unique(All(ind).out.exp.stimID);
    
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    v = 1; % no vis stim condition...
    
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    
    for i = 1:numStims
        
        holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
        c=c+1;
       
        
        mRespEns{c} = squeeze(respMat(i, v, :) - ...
           baseMat(i, v, :));
        
        if i==1
            mRespEns{c}(ROIinArtifact) = nan;
        else
%             c2 = c2+1;
%             if c2==659
%                 disp('Hi there')
%             end
            mRespEns{c}(offTargetRisk(holo, :)' | ROIinArtifact) = nan; % was 'and' but bc we're removing here, should not be an or asi is elsewhere in code
        end
        
        isControl(c) = i==1; % no stim condition
        uExpts(c) = ind;
    end 
end

mRespEns = mRespEns(~isControl);
uExpts = uExpts(~isControl);

outVars.mRespEns = mRespEns;
outVars.uExptsEns = uExpts;