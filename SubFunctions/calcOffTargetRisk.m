function [All] = calcOffTargetRisk(All,opts)

thisPlaneTolerance = opts.thisPlaneTolerance;
onePlaneTolerance = opts.onePlaneTolerance;
numExps = numel(All);

for ind =1:numExps
    %%offtargetRisk
    stimCoM = All(ind).out.exp.stimCoM;
    numCells = size(All(ind).out.exp.allData,1);
    allCoM = All(ind).out.exp.allCoM;
    stimDepth = All(ind).out.exp.stimDepth;
    allDepth = All(ind).out.exp.allDepth;
    muPerPx = 800/512;
    
    allLoc = [allCoM*muPerPx (allDepth-1)*30];
    stimLoc = [stimCoM*muPerPx (stimDepth-1)*30];
    
    roisTargets = All(ind).out.exp.rois;
    holoTargets = All(ind).out.exp.holoTargets;
    
    radialDistToStim=zeros([size(stimCoM,1) numCells]);
    axialDistToStim = zeros([size(stimCoM,1) numCells]);
    StimDistance = zeros([size(stimCoM,1) numCells]);
    for i=1:size(stimCoM,1)
        for k=1:numCells
            D = sqrt(sum((stimCoM(i,:)-allCoM(k,:)).^2));
            radialDistToStim(i,k)=D;
            z = stimDepth(i)-allDepth(k);
            axialDistToStim(i,k) = z;
            StimDistance(i,k) = sqrt(sum((stimLoc(i,:)-allLoc(k,:)).^2));
        end
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
            temp3 = radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==2;
        else
            temp = any(radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0);
            temp2 = any(radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1);
            temp3 = any(radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==2);
        end
        offTargetRisk(i,:) = temp | temp2 | temp3;
    end
    All(ind).out.anal.offTargetRisk = offTargetRisk;
end

end

