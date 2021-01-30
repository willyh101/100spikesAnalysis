function All = calcDeltaActivityChettihHarveyBAD(All)
% influence measurement
% Chettih and Harvey (2019) "Single-neuron perturbations reveal
% feature-specific competition in V1"
% https://www.nature.com/articles/s41586-019-0997-6
%
% technically, this is the deconvolved implentation of thier algorithm, it
% was also used on dF data, but for dF they used a vector of time points
% aligned to stimulus onset instead of the single scalar value of
% single-trial activity.


for ind = 1:numel(All)
    v = 1; % restrict to no vis for now
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;

    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);
    ctrl = us(1);
    
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;

    % mean control response is define as the mean across all no holo and no
    % vis stim trials (baseline subtracted.... rdData - bdata)
    mean_control_response = mean(All(ind).out.exp.rdData(:,...
        trialsToUse & All(ind).out.exp.stimID == ctrl &...
        All(ind).out.exp.visID == v), 2)...      
        - mean(All(ind).out.exp.bdata(:, trialsToUse & ...
        All(ind).out.exp.stimID == ctrl &...
        All(ind).out.exp.visID == v), 2);
    
    clear mActivity dActivity sActivity rActivity
    
    
    for i=1:numel(us)
        s = us(i);
        
        if i==1
            cellsToUse = ~ROIinArtifact';
        else
            cellsToUse = ~ROIinArtifact' & ~offTargetRisk(i-1,:);
        end
        
        % get the baseline subtracted response
        data = All(ind).out.exp.rdData - All(ind).out.exp.bdata;
        
        % mActivity -> the mean of a cell x trial's response - the mean
        % control response for a cell.
        mActivity(:,i) = mean(data(:,trialsToUse ...
             & All(ind).out.exp.stimID == s)...
            - mean_control_response, 2); 
        
        % sActivity -> STD of those responses
        sActivity(:,i) = std(data(:,trialsToUse ...
             & All(ind).out.exp.stimID == s)...
            - mean_control_response, [], 2);
        
        % dActivity -> "delta activity", the CH influence metric
        dActivity(:,i) = mActivity(:,i)./sActivity(:,i);
        
        % rActivity -> single trial residuals
        rActivity{i} = data(:,trialsToUse ...
             & All(ind).out.exp.stimID == s)...
            - mean_control_response;
        
        % TESTING
        mean_s = rActivity{i}./std(rActivity{i},[],2);
        mean_r = mean(mean_s,2);
        
        mActivity(~cellsToUse,i) = nan;
        sActivity(~cellsToUse,i) = nan;
        dActivity(~cellsToUse,i) = nan;
        rActivity{i}(~cellsToUse,:) = nan;
        
    end
    
    All(ind).out.anal.dActivity = dActivity; % delta-activity the main output (mean/std)
    All(ind).out.anal.sActivity = sActivity; % std of difference in activity across all trials
    All(ind).out.anal.mActivity = mActivity; % mean difference of stims by trial 
    All(ind).out.anal.rActivity = rActivity; % single trial residuals (the un-meaned mActivity), a cell array bc there are diff num trials for each stim
    All(ind).out.anal.cActivity = mean_control_response; % the mean control, no response condition
    
end