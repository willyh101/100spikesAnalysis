%%
% New pre-processing step (added May, 2022)
%
% Remove cells from analysis on the current trial if they were an 
% offTarget on the previous trial 
%%
function [All] = prevOffTargetRemoval(All)

numExps = numel(All);
%% Loop through all experiments
for ind=1:numExps
    unIDs = unique(All(ind).out.exp.stimID);
    totalCells = size(All(ind).out.exp.dfData,1);
    totalTrials = size(All(ind).out.exp.stimID,2);
     
    % Everyone starts being avaliable for analysis
    bufferZone = 1;
    lastTimeStimmed = bufferZone+zeros(totalCells,1);
    % Loop through all trials sequentially
    for ii = 1:totalTrials
        % Remove these cells from this trial
        doNotInclude = lastTimeStimmed<bufferZone;
        
        % Update the lastTimeStimmed list
        % Figure out who is currently stimmed
        tempIndex = find(unIDs==All(ind).out.exp.stimID(ii),1);
        holoIndex = All(ind).out.exp.stimParams.roi{tempIndex};
        if holoIndex > 0
            currTargets = All(ind).out.exp.holoTargets{holoIndex}(~isnan(All(ind).out.exp.holoTargets{holoIndex}));
            currOffTargets=find(All(ind).out.anal.offTargetRisk(holoIndex,:)==1);
        else
            % no stim, so no current targets
            currTargets = [];
            currOffTargets=[];
        end
        lastTimeStimmed = lastTimeStimmed + 1;
        % Reset the count of those just stimmed
%         lastTimeStimmed(currOffTargets) = 0;
        lastTimeStimmed(currTargets) = 0;
        
        % Cells that are currently being stimmed are not removed 
        % (needed to determine if the stim was successful later on)
        doNotInclude(currTargets) = 0;
        
        % Eliminate these cells/trials from data
        All(ind).out.exp.dfData(doNotInclude,:,ii) = nan;
        All(ind).out.exp.zdfData(doNotInclude,:,ii) = nan;
        All(ind).out.exp.allData(doNotInclude,:,ii) = nan;
    end
    
    
    %% Check to make sure each cell is still in at least 10 trials 
    % (and remove those that don't)
    
    % Loop through all of the stims
    for gh=1:length(unIDs)
        
        tempIndex = gh;
        holoIndex = All(ind).out.exp.stimParams.roi{tempIndex};
        
        trialsToUse = ismember(All(ind).out.exp.stimID, unIDs(gh)) &...
            All(ind).out.exp.lowMotionTrials &...
            All(ind).out.exp.lowRunTrials & ...
            All(ind).out.exp.stimSuccessTrial & ...
            All(ind).out.exp.visID==1;
                
        tempTotalTrials = sum(trialsToUse);
        badCells = (tempTotalTrials-sum(isnan(All(ind).out.exp.dfData(:,1,trialsToUse)),3))<10;
       
        % Cells that are currently being stimmed are not removed
        if holoIndex > 0
            badCells(All(ind).out.exp.holoTargets{holoIndex}(~isnan(All(ind).out.exp.holoTargets{holoIndex})))=0;
        end
                
        All(ind).out.exp.dfData(badCells,:,trialsToUse) = nan;
        All(ind).out.exp.zdfData(badCells,:,trialsToUse) = nan;
        All(ind).out.exp.allData(badCells,:,trialsToUse) = nan;             
    end
end

end

