function [outVars] = grandDistanceMaker(opts,All,outVars)
distType = opts.distType;
numExps = numel(All);

distToEnsemble=[];
offTargetRiskEns=[];
for ind = 1:numExps
    
    %First create stimDistance
    if ~isfield(All(ind).out.anal,'StimDistance')
        
        stimCoM = All(ind).out.exp.stimCoM;
        numCells = size(All(ind).out.exp.zdfData,1);
        allCoM = All(ind).out.exp.allCoM;
        stimDepth = All(ind).out.exp.stimDepth;
        allDepth = All(ind).out.exp.allDepth;
        muPerPx = 800/512;
        
        allLoc = [allCoM*muPerPx (allDepth-1)*30];
        stimLoc = [stimCoM*muPerPx (stimDepth-1)*30];
        
        StimDistance = zeros([size(stimCoM,1) numCells]);
        for i=1:size(stimCoM,1);
            for k=1:numCells;
                StimDistance(i,k) = sqrt(sum((stimLoc(i,:)-allLoc(k,:)).^2));
            end
        end
    else
        StimDistance = All(ind).out.anal.StimDistance;
    end
    
    %load in some variables you'll need
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    
    %     now iterate through every stim and see the response as a function of
    %     distance
    us = unique(All(ind).out.exp.stimID);
    for k=1:numel(us)
        s = us(k);
        sidx = find(us==s);
        holo = All(ind).out.exp.stimParams.roi{sidx}; % Better Identifying ensemble
        
        if holo>0
            
            Tg=All(ind).out.exp.rois{holo};
            dists = StimDistance(Tg,:);
                   
            switch distType
                case 'min'
                    minDist = min(dists,[],1);
                    distToUse = minDist;
                case 'max'
                    maxDist = max(dists,[],1);
                    distToUse = maxDist;
                case 'geo'
                    try
                        geoDist = geo_mean(dists,1);
                    catch
                        geoDist = geomean(dists,1);
                    end
                    distToUse = geoDist;
                case 'mean'
                    meanDist = mean(dists,1);
                    distToUse = meanDist;
                case 'harm'
                    harmDist = harmmean(dists,1);
                    distToUse = harmDist;
                otherwise
                    disp('dist type not understood, using min')
                    minDist = min(dists,[],1);
                    distToUse =minDist;
            end
            
            distToEnsemble{end+1} = distToUse;
            offTargetRiskEns{end+1} =  ROIinArtifact' |...
                offTargetRisk(holo,:);
        end
    end
end
outVars.distToEnsemble = distToEnsemble;
outVars.offTargetRiskEns = offTargetRiskEns;