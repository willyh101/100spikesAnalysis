function [popRespDistEns] = popDistMaker_GH(opts,All,CellToUseVar,plotFlag)
%% make distance plots 
distType = lower(opts.distType);
distBins = opts.distBins; %[0:25:1000];
numDist =numel(distBins)-1;
numExps = numel(All);

c=0;
for ind = 1:numExps
    %add additional constraints to cellsToUse
    if isempty(CellToUseVar)
        cellToUseLimit = ones([1 size(All(ind).out.anal.respMat,3)]);
    elseif islogical(CellToUseVar)
        cellToUseLimit=CellToUseVar;
    else
       try
            cellToUseLimit = eval(['All(ind).out.' CellToUseVar]);
       catch
            disp('ERROR variable not found')
            break
        end
    end
    %First create stimDistance
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
    
    %load in some variables you'll need
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    
%     now iterate through every stim and see the response as a function of
%     distance
    for i= 1:numStims
        holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
        
        if i~=1
            c=c+1;
            Tg=All(ind).out.exp.rois{holo};
            dists = StimDistance(Tg,:);
            
            
            center_of_mass = mean(stimLoc(All(ind).out.exp.rois{holo},:));
%             centerDist = sqrt(sum((center_of_mass-allLoc).^2,2))';
%             centerDistbyHolo(i,:) = centerDist;
                        
%             minDistbyHolo(i,:) = minDist;
%             geoDistbyHolo(i,:) = geoDist;
%             meanDistbyHolo(i,:) = meanDist;
%             harmDistbyHolo(i,:) = harmDist;
            
            switch distType
                case 'min'
                    minDist = min(dists,[],1);
                    distToUse = minDist;
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
                case 'center of mass'
                    centerDist = sqrt(sum((center_of_mass-allLoc).^2,2))';
                    distToUse = centerDist;
                otherwise
                    disp('dist type not understood, using min')
                    minDist = min(dists,[],1);
                    distToUse =minDist;
            end
            
            
            
            for d = 1:numel(distBins)-1
                cellsToUse = ~ROIinArtifact' &...
                    ~offTargetRisk(holo,:) &...
                    distToUse > distBins(d) &...
                    distToUse <= distBins(d+1) &...
                    cellToUseLimit;
                %if you eventually want to use this with different vis
                %stims use the commented out code.
%                 popRespDistEns(c,v,d) = nanmean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                popRespDistEns(c,d) = nanmean(squeeze(respMat(i,1,cellsToUse) - baseMat(i,1,cellsToUse)));
            end
        end
    end
    
%     popRespDistEns{ind}=popRespDist;
end
if plotFlag
    figure(4);clf
    imagesc(popRespDistEns)
    colormap rdbu
end
caxis([-0.3 +0.3])