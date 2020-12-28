function All = shuffleChettihHarvey(All)

nboot = 100000;

for ind = 1:numel(All)

    v = 1;
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;

    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);
    ctrl = us(1);

    shuffledStims = All(ind).out.exp.stimID(randperm(length(All(ind).out.exp.stimID)));


    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;


    mean_control_response = mean(All(ind).out.exp.rdData(:,...
        trialsToUse & shuffledStims == ctrl &...
        All(ind).out.exp.visID == v), 2)...      
        - mean(All(ind).out.exp.bdata(:, trialsToUse & ...
        shuffledStims == ctrl &...
        All(ind).out.exp.visID == v), 2);

    clear mActivity dActivity sActivity rActivity 
        
        
        for i=1:numel(us)
            activityBoot = zeros(length(ROIinArtifact),nboot);
            parfor b=1:nboot
                disp(b)
                if mod(b,100) == 0
                    fprintf('.')
                end
                if mod(b,10000) == 0
                    fprintf('\r')
                end
                s = us(i);

                shuffledStims = All(ind).out.exp.stimID(randperm(length(All(ind).out.exp.stimID)));

                if i==1
                    cellsToUse = ~ROIinArtifact';
                else
                    % do any here bc the stims aren't matched to targets
                    cellsToUse = ~ROIinArtifact' & ~any(offTargetRisk); 
                end


                data = All(ind).out.exp.rdData - All(ind).out.exp.bdata;
                mActivity = mean(data(:,trialsToUse ...
                     & shuffledStims == s)...
                    - mean_control_response, 2); % aka single trial residuals
                sActivity = std(data(:,trialsToUse ...
                     & All(ind).out.exp.stimID == s)...
                    - mean_control_response, [], 2);
                dActivity = mActivity(:,i)./sActivity(:,i);


                dActivity(~cellsToUse,:) = nan;


                activityBoot(:,b) = dActivity;
            end
            activityBoot = mean(activityBoot,1);
            dActivityshuf(:,i) = activityBoot;
        end
end

All(ind).out.anal.dActivityShuf = dActivityshuf; % delta-activity the main output (mean/std)

    
end