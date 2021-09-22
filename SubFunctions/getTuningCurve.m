function [All, outVars] = getTuningCurve(All, opts, outVars)

% general settings
visAlpha = opts.visAlpha;
numExps = numel(All);

for ind=1:numExps  
    uVisID = unique(All(ind).out.vis.visID);
    uVisID(uVisID==0)=[];
    
    oriCurve=[];
    oriCurveSEM = [];
    oriTrialCount =[];
    for i=1:numel(uVisID)
        v= uVisID(i);
        
        trialsToUse = All(ind).out.vis.visID==v &...
            All(ind).out.vis.lowMotionTrials &...
            All(ind).out.vis.lowRunTrials;
        
        oriCurve(i,:) = mean(All(ind).out.vis.rdData(:,trialsToUse), 2);
        oriCurveSEM(i,:) = sem2(All(ind).out.vis.rdData(:,trialsToUse), 2);
        oriTrialCount(i) = sum(trialsToUse); 
    end
    
    All(ind).out.anal.oriCurve = oriCurve;
    All(ind).out.anal.oriCurveSEM = oriCurveSEM;
    
    [~, maxOriIndex]= max(oriCurve);
    All(ind).out.anal.prefOri = maxOriIndex;
    
    prefOri = maxOriIndex;
    orthoOri = prefOri-2;
    orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
    
    orthoOri2 = orthoOri+4;
    orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
    
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
    