function [popRespCorr] = popCorrMakerSingle(opts,All,CellToUseVar,plotFlag,ensID)



distBins = opts.CorrSpace;
stringCorrType = opts.CorrToPlot;

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

vs =  unique(All(ind).out.exp.visID);
vs(vs==0)=[];
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
    Tg=All(ind).out.exp.holoTargets{holo};
    Tg(isnan(Tg))=[];
    
    distCorr = corrToUse(Tg,:);
    if numel(Tg)==1
        minDist = distCorr;
    else
        minDist = nanmean(distCorr);
    end
    
    if numel(Tg)==0
        minDistbyHolo(i,:) = ones([1 size(minDist,2)])*1000;
    else
        minDistbyHolo(i,:) = minDist;
    end
    %                 distBins = linspace(-0.5,0.5,40);
    for d = 1:numel(distBins)-1
        cellsToUse = ~ROIinArtifact' &...
            ~offTargetRisk(holo,:) &...
            minDist > distBins(d) &...
            minDist <= distBins(d+1) &...
            cellToUseLimit;
        
        popRespCorr(d) = nanmean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
        
        noHoloEquivalent = nanmean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse)));
        popRespCorrSub(d) =  popRespCorr(d) - noHoloEquivalent;
    end
else
    
    popRespCorr =nan([numel(distBins)-1 1]);
end

if plotFlag
    figure(5);clf
    p =plot(popRespCorr);
    p.LineWidth=2;
    ylim([-0.3 +0.3])
end
