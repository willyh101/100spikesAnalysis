%%offtargetRisk
thisPlaneTolerance = opts.thisPlaneTolerance;
onePlaneTolerance = opts.onePlaneTolerance;

offTargetRiskEns=[];

disp('Recalculating OffTargetRisk')
for ind =1:numExps
    fprintf('.')
    dataToUse = All(ind).out.exp.dataToUse;
    
    stimCoM = All(ind).out.exp.stimCoM;
    numCells = size(dataToUse,1);
    allCoM = All(ind).out.exp.allCoM;
    stimDepth = All(ind).out.exp.stimDepth;
    allDepth = All(ind).out.exp.allDepth;
    muPerPx = 800/512;
      
    roisTargets = All(ind).out.exp.rois;
    holoTargets = All(ind).out.exp.holoTargets;
    
    if ~isfield(All(ind).out.anal, 'radialDistToStim')
        radialDistToStim=zeros([size(stimCoM,1) numCells]);
        axialDistToStim = zeros([size(stimCoM,1) numCells]);
        for i=1:size(stimCoM,1);
            for k=1:numCells;
                D = sqrt(sum((stimCoM(i,:)-allCoM(k,:)).^2));
                radialDistToStim(i,k)=D;
                z = stimDepth(i)-allDepth(k);
                axialDistToStim(i,k) = z;
            end
        end
        All(ind).out.anal.radialDistToStim = radialDistToStim;
        All(ind).out.anal.axialDistToStim = axialDistToStim;
    else
        radialDistToStim = All(ind).out.anal.radialDistToStim;
        axialDistToStim = All(ind).out.anal.axialDistToStim;
    end
    
    offTargetRisk = zeros([numel(roisTargets) numCells]);
    for i=1:numel(roisTargets)
        Tg = roisTargets{i};
        try
            TgCells = holoTargets{i};
        catch;end;
        
        if numel(Tg) == 1
            temp = radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0;
            temp2 = radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1;
        else
            temp = any(radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0);
            temp2 = any(radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1);
        end
        offTargetRisk(i,:) = temp | temp2;
    end
    All(ind).out.anal.offTargetRisk = offTargetRisk;
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    us = unique(All(ind).out.exp.stimID);
    for k=1:numel(us)
        s = us(k);
        sidx = find(us==s);
        holo = All(ind).out.exp.stimParams.roi{sidx}; % Better Identifying ensemble
        if holo~=0
        offTargetRiskEns{end+1} =  ROIinArtifact' |...
            offTargetRisk(holo,:);
        end
    end
end
outVars.offTargetRiskEns=offTargetRiskEns;


disp('done')