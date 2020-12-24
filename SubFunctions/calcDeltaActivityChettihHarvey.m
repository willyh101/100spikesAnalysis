function All = calcDeltaActivityChettihHarvey(All)

for ind = 1:numel(All)
    v = 1;
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;

    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);
    ctrl = us(1);
    
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;

    
    mean_control_response = mean(All(ind).out.exp.rdData(:,...
        trialsToUse & All(ind).out.exp.stimID == ctrl &...
        All(ind).out.exp.visID == v), 2)...      
        - mean(All(ind).out.exp.bdata(:, trialsToUse & ...
        All(ind).out.exp.stimID == ctrl &...
        All(ind).out.exp.visID == v), 2);
    
    clear mActivity dActivity sActivity
    for i=1:numel(us)
        s = us(i);
        
        if i==1
            cellsToUse = ~ROIinArtifact';
        else
            cellsToUse = ~ROIinArtifact' & ~offTargetRisk(i-1,:);
        end
        
        
        data = All(ind).out.exp.rdData - All(ind).out.exp.bdata;
        mActivity(:,i) = mean(data(:,trialsToUse ...
             & All(ind).out.exp.stimID == s)...
            - mean_control_response, 2); % aka single trial residuals
        sActivity(:,i) = std(data(:,trialsToUse ...
             & All(ind).out.exp.stimID == s)...
            - mean_control_response, [], 2);
        dActivity(:,i) = mActivity(:,i)./sActivity(:,i);
        
        mActivity(~cellsToUse,i) = nan;
        sActivity(~cellsToUse,i) = nan;
        dActivity(~cellsToUse,i) = nan;
        
    end
    
    All(ind).out.anal.dActivity = dActivity; % delta-activity the main output (mean/std)
    All(ind).out.anal.sActivity = sActivity; % std of difference in activity across all trials
    All(ind).out.anal.mActivity = mActivity; % mean difference of stims by trial
    All(ind).out.anal.cActivity = mean_control_response; % the mean control, no response condition
    
end