function All = calcDeltaActivityChettihHarvey(All)
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
    
    clear mInfluence rInfluence
    
    for i=1:numel(us)
        s = us(i);
        
        if i==1
            cellsToUse = ~ROIinArtifact';
        else
            cellsToUse = ~ROIinArtifact' & ~offTargetRisk(i-1,:);
        end
        
        % baseline data first
        data = All(ind).out.exp.rdData - All(ind).out.exp.bdata;
        
        % calculate single trial response residuals (to control)
        resid = data(:, trialsToUse & All(ind).out.exp.stimID == s)...
            - mean_control_response;
        
        % normalized cell-by-cell to stf of those residuals
        normalized_residuals = resid./std(resid,[],2);
        
        % remove off targets and artifact cells
        normalized_residuals(~cellsToUse,:) = nan;
        
        % mean influence is the cellwise mean
        mInfluence(:,i) = mean(normalized_residuals,2);
        rInfluence{i} = normalized_residuals(cellsToUse,:);
        
    end
    
    All(ind).out.anal.mInfluence = mInfluence;
    All(ind).out.anal.rInfluence = rInfluence;
end