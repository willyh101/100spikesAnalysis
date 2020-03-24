function [All, outVars] = getTuningCurve(All, opts, outVars)

% general settings
visAlpha = opts.visAlpha;
numExps = numel(All);

for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    uVisID = unique(All(ind).out.vis.visID);
    uVisID(uVisID==0)=[];
    
    oriCurve=[];
    for i=1:numel(uVisID)
        v= uVisID(i);
        
        trialsToUse = All(ind).out.vis.visID==v &...
            All(ind).out.vis.lowMotionTrials &...
            All(ind).out.vis.lowRunTrials;
        
        oriCurve(i,:)=mean(All(ind).out.vis.rdata(:,trialsToUse),2);
    end
    
    All(ind).out.anal.oriCurve = oriCurve;
    
    [~, maxOriIndex]= max(oriCurve);
    All(ind).out.anal.prefOri = maxOriIndex;
    
    prefOri = maxOriIndex;
    orthoOri = prefOri-2;
    orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
    
    orthoOri2 = orthoOri+4;
    orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
    
    orthoOri = cat(1,orthoOri, orthoOri2);
    
    % save the datas
    prefOris{ind} = prefOri;
    orthoOris{ind} = orthoOri;
    All(ind).out.anal.prefOri = prefOri;
    All(ind).out.anal.orthoOri = orthoOri;

end

outVars.prefOris = prefOris;
outVars.orthoOris = orthoOris;
    