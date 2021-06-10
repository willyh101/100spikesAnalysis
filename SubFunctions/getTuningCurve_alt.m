function [All, outVars] = getTuningCurve_alt(All, opts, outVars)

% general settings
visAlpha = opts.visAlpha;
numExps = numel(All);
os = [nan 0:45:135];
for ind=1:numExps  
    uVisID = [1 2 3 4 5];
    
    oriCurve=[];
    oriCurveSEM = [];
    oriTrialCount =[];
    
    for i=1:numel(uVisID)
        v= uVisID(i);
        
        visIDs = All(ind).out.vis.visID;
        visIDs(visIDs>5) = visIDs(visIDs>5)-4;
        
        trialsToUse = visIDs==v &...
            All(ind).out.vis.lowMotionTrials &...
            All(ind).out.vis.lowRunTrials;
        
        oriCurve(i,:) = mean(All(ind).out.vis.rdata(:,trialsToUse), 2);
        oriCurveSEM(i,:) = sem2(All(ind).out.vis.rdata(:,trialsToUse), 2);
        oriTrialCount(i) = sum(trialsToUse); 
    end
    
    All(ind).out.anal.oriCurve = oriCurve;
    All(ind).out.anal.oriCurveSEM = oriCurveSEM;
    
    [~, maxOriIndex]= max(oriCurve);
    All(ind).out.anal.prefOri = maxOriIndex;
    
    prefOri = maxOriIndex;
    orthoOri = prefOri-2;
    orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+4;
    
    orthoOri2 = orthoOri+4;
    orthoOri2(orthoOri2>5) = orthoOri2(orthoOri2>5)-4;
    
    orthoOri = cat(1,orthoOri, orthoOri2);
    
    % save the datas
    tuningCurves{ind} = oriCurve;
    tuningCurvesSEM{ind} = oriCurveSEM;
    tuningCurveTrials{ind} = oriTrialCount;
    prefOris{ind} = prefOri;
    orthoOris{ind} = orthoOri;
    isVisR{ind} = All(ind).out.anal.pVisR < 0.05;

    
    All(ind).out.anal.prefOri = prefOri;
    All(ind).out.anal.orthoOri = orthoOri;

end

outVars.prefOris = prefOris;
outVars.orthoOris = orthoOris;
outVars.tuningCurves = tuningCurves;
outVars.tuningCurvesSEM = tuningCurvesSEM;
outVars.tuningCurveTrials = tuningCurveTrials;
outVars.isVisR = isVisR;
    