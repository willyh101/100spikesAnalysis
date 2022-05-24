%%
% Estimates whether the stim successfully evoked a response in the 
% targeted cells
%%
function [All] = stimSuccessFn(All, opts)

numExps = numel(All);

stimsuccessZ = opts.stimsuccessZ;%0.25; %over this number is a succesfull stim
stimEnsSuccess = opts.stimEnsSuccess;% 0.5; %fraction of ensemble that needs to be succsfull

for ind =1:numExps
    
    
    opts.recWinRange = [0.5 1.5];
    recWinRange =opts.recWinRange;%
    visStart = find(diff(All(ind).out.exp.outputsInfo.OutputPatterns{1}(:,9)>0),1)/20000;
    if isempty(visStart)
        fprintf('\nNo Exp VisStart Detected...');
        visStart=0.5;
    end
    
    sz = size(All(ind).out.exp.dfData);
    
    recWinSec = recWinRange + visStart;
    
    if ~isfield(All(ind).out.info,'FR')
        All(ind).out.info.FR=6;
    end
    
    winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
    bwinToUse = max(floor([0 visStart]*All(ind).out.info.FR),[1 1]);
        
    %ensure has a visID
    if ~isfield(All(ind).out.exp,'visID')
        All(ind).out.exp.visID = ones(size(All(ind).out.exp.stimID));
        disp(['Added visID to Exp ' num2str(ind)]);
    end
    
    % create a trial by trial stimSuccess limit
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.visID);
%     numTrials = size(dataToUse,3);
    numTrials = size(All(ind).out.exp.stimID,2);
    
    if isfield(opts, 'stimSuccessByZ') && opts.stimSuccessByZ
        rdata = squeeze(nanmean(All(ind).out.exp.zdfData(:,winToUse(1):winToUse(2),:),2));
        bdata = squeeze(nanmean(All(ind).out.exp.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
    else
        rdata = All(ind).out.exp.rdData;
        bdata = All(ind).out.exp.bdata;
    end
    
    %ensure stimparam correct and properly formatted
    %Caution may erase stimparams if they are complex
    for r = 1:size(All(ind).out.exp.stimParams.roi,2)
        x = All(ind).out.exp.stimParams.roi{r};
        if iscell(x)
            All(ind).out.exp.stimParams.roi{r} = x{1};
        end
        if numel(x)>1
            All(ind).out.exp.stimParams.roi{r} = r-1;
        end
    end
    
    clear stimSuccessTrial
    for k=1:numTrials
        s = All(ind).out.exp.stimID(k);
        sidx = find(us==s);
        v = All(ind).out.exp.visID(k);
        vidx = find(vs==v);
        
        h = All(ind).out.exp.stimParams.roi{sidx};
        if h==0
            htg=[];
            stimSuccessTrial(k) = 1;
        else
            try
                htg = All(ind).out.exp.holoTargets{h};
            catch
                disp('Incorrect number of holoTargets')
                htg =[];
            end
            htg(isnan(htg))=[];
            
            vals = rdata(htg,k) - bdata(htg,k);
            stimScore = vals>stimsuccessZ;
            stimSuccessTrial(k)= mean(stimScore) > stimEnsSuccess;
            tempStimScore(k) = mean(stimScore);
        end
    end
    
    
    All(ind).out.exp.stimSuccessTrial = stimSuccessTrial;
    percentSuccessStim(ind)=mean(stimSuccessTrial);
    
    
end


end

