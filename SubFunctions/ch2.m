for ind = 1%:numel(All)
    v = 1; % restrict to no vis for now
    visToUse = All(ind).out.exp.visID == 1 | All(ind).out.exp.visID == 0;
    
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        visToUse;

    us = unique(All(ind).out.exp.stimID);
    ctrl = All(ind).out.exp.stimID == us(1);
    
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    cellsToUse = ~ROIinArtifact' & ~any(offTargetRisk);
    

    % mean control response is define as the mean across all no holo and no
    % vis stim trials (baseline subtracted.... rdData - bdata)
    data = All(ind).out.exp.rdData- All(ind).out.exp.bdata;
    
    mean_control_response = mean(data(cellsToUse, trialsToUse & ctrl), 2);

    all_other_responses = data(cellsToUse, trialsToUse & ~ctrl);
    
    single_trial_residuals = all_other_responses - mean_control_response;
    
    normalized_residuals = single_trial_residuals./std(single_trial_residuals,[],2);
    
    for i=1:numel(us)
        s = us(i);
        stims = All(ind).out.exp.stimID == s;
    
    
end
    
    