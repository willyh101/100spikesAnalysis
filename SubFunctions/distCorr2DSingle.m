function [popRespCorr numResponders] = distCorr2DSingle(opts,All,CellToUseVar,plotFlag,ensID)



distBinsCorr = opts.CorrSpace;
stringCorrType = opts.CorrToPlot;

distType = opts.distType;
distBins = opts.distBins; %[0:25:1000];
numDist =numel(distBins)-1;
numExps = numel(All);

%% Plot
% disp('Calculating...')

clear popResponseCorr

ind=1;

if isempty(CellToUseVar)
    cellToUseLimit = ones([1 size(All(ind).out.anal.respMat,3)]);
elseif islogical(CellToUseVar)
    cellToUseLimit=CellToUseVar;
else
    try
        cellToUseLimit = eval(['All(ind).out.' CellToUseVar]);
    catch
        disp('ERROR CellToUseVar not found')
        return
    end
end

try
    corrToUse  = eval(['All(ind).out.anal.' stringCorrType]); %AllCorr;
catch
    disp('That Corr type not available')
    return
end

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

respMat = All(ind).out.anal.respMat;
baseMat = All(ind).out.anal.baseMat;

ROIinArtifact = All(ind).out.anal.ROIinArtifact;
offTargetRisk = All(ind).out.anal.offTargetRisk;

clear popRespCorr minDistbyHolo cellsToUse
v=1;
i=ensID;
holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
%             popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));

if i~=1
    Tgc=All(ind).out.exp.holoTargets{holo};
    Tgc(isnan(Tgc))=[];
    
    Tgd=All(ind).out.exp.rois{holo};
    dists = StimDistance(Tgd,:);

    
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
        otherwise
            disp('dist type not understood, using min')
            minDist = min(dists,[],1);
            distToUse =minDist;
    end
    
    
    distCorr = corrToUse(Tgc,:);
    if numel(Tgc)==1
        minDistCorr = distCorr;
    else
        minDistCorr = nanmean(distCorr);
    end
    
    
    %                 distBins = linspace(-0.5,0.5,40);
    for dc = 1:numel(distBinsCorr)-1
        for dd = 1:numel(distBins)-1
            cellsToUse = ~ROIinArtifact' &...
                ~offTargetRisk(holo,:) &...
                minDistCorr > distBinsCorr(dc) &...
                minDistCorr <= distBinsCorr(dc+1) &...
                distToUse > distBins(dd) &...
                distToUse <= distBins(dd+1) &...
                cellToUseLimit;
            
            popRespCorr(dc,dd) = nanmean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            noHoloEquivalent = nanmean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse)));
            popRespSub(dc,dd) =  popRespCorr(dc) - noHoloEquivalent;
            numResponders(dc,dd) = sum(cellsToUse);
        end
    end
    
else
    
    popRespCorr = nan([numel(distBinsCorr)-1 numel(distBins)-1]);
    numResponders = nan([numel(distBinsCorr)-1 numel(distBins)-1]);
end

if plotFlag
    figure(5);clf
    p =plot(popRespCorr);
    p.LineWidth=2;
    ylim([-0.3 +0.3])
end
