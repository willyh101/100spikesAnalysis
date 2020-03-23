function [All,outVars] = cleanData(All,opts)
%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=opts.FRDefault;%6;
recWinRange =opts.recWinRange;% [0.5 1.5];% %from vis Start [1.25 2.5];


%Stim Success Thresholds
stimsuccessZ = opts.stimsuccessZ;%0.25; %over this number is a succesfull stim
stimEnsSuccess = opts.stimEnsSuccess;% 0.5; %fraction of ensemble that needs to be succsfull

%run Threshold
runThreshold = opts.runThreshold;% 6 ; %trials with runspeed below this will be excluded

numExps = numel(All);

clear ensStimScore

 for ind =1:numExps
     pTime =tic;
     fprintf(['Processing Experiment ' num2str(ind) '...']);
     
     All(ind).out.anal.numCells = size(All(ind).out.exp.zdfData,1);
     numCells(ind) = size(All(ind).out.exp.zdfData,1);
     
     if ~isfield(All(ind).out.info,'FR')
         All(ind).out.info.FR=FRDefault;
     end
     
     sz = size(All(ind).out.exp.zdfData);
         
      try
        visStart = find(diff(All(ind).out.exp.outputsInfo.OutputPatterns{1}(:,9)>0),1)/20000;
    catch
        fprintf('\nError in detecting VisStart...');
        visStart = 0.5;
    end
    if isempty(visStart);
        fprintf('\nNo Exp VisStart Detected...');
        visStart=0.5;
    end
     All(ind).out.exp.visStart = visStart;
%     
%      if ~isfield(All(ind).out.exp,'visStart')
%          All(ind).out.exp.visStart = 0.5;
%          disp(['Added visStart to Exp' num2str(ind)]);
%      end
%      
     recWinSec = recWinRange + All(ind).out.exp.visStart; %recording window relative to when vis start
     
     
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
     bwinToUse = max(round([0 All(ind).out.exp.visStart]*All(ind).out.info.FR),[1 1]);

     rdata = squeeze(mean(All(ind).out.exp.zdfData(:,winToUse(1):winToUse(2),:),2));
     bdata = squeeze(mean(All(ind).out.exp.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
     
     All(ind).out.exp.rdData=rdata;
     All(ind).out.exp.bdata=bdata;
     
     % hacky fix for exp11 with no vis stim
     if ind==11
         vrdata = squeeze(mean(All(ind).out.vis.zdfData(:,winToUse(1):winToUse(2),:),2));
         vbdata = squeeze(mean(All(ind).out.vis.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
         All(ind).out.vis.rdData=vrdata;
         All(ind).out.vis.bdata=vbdata;
     end
         
     
     All(ind).out.anal.recWinUsed = winToUse;
    All(ind).out.anal.bwinToUse = bwinToUse;
    All(ind).out.anal.recStartTime = visStart;
    All(ind).out.anal.recStartFrame = round(visStart*All(ind).out.info.FR);
  
     %runProcessing Section
    runVal = All(ind).out.exp.runVal;
    rnSz = size(runVal);
    runperiod = [1:min(All(ind).out.anal.recWinUsed(2),rnSz(2))];
   
    lowRunVals = mean((runVal(:,runperiod)<runThreshold)');
    lowRunTrials = lowRunVals>0.75; %percent of frames that need to be below run threshold
    if numel(lowRunTrials)>numel(All(ind).out.exp.stimID)
        lowRunTrials = lowRunTrials(1:numel(All(ind).out.exp.stimID));
    end
    All(ind).out.exp.lowRunTrials = lowRunTrials;
    
    percentLowRunTrials(ind) = mean(lowRunTrials);

    
    if ~isfield(All(ind).out.vis,'lowRunTrials')
        try
        runVal = All(ind).out.vis.runVal;
        rnSz = size(runVal);
        runperiod = [1:min(All(ind).out.anal.recWinUsed(2),rnSz(2))];
        lowRunVals = mean((runVal(:,runperiod)<runThreshold)');
        lowRunTrials = lowRunVals>0.75; %percent of frames that need to be below run threshold
        All(ind).out.vis.lowRunTrials = lowRunTrials;
        catch
            disp('Vis Run Trial Problem')
        end
    end
       
    
    %Total Number of Targets Shot per recording
     temp = unique([All(ind).out.exp.holoTargets{:}]);
     temp(isnan(temp))=[];
     All(ind).out.anal.targets = temp;
     numUniqueTargets(ind) =numel(temp);
     
     %ensure has a visID
     if ~isfield(All(ind).out.exp,'visID')
         All(ind).out.exp.visID = ones(size(All(ind).out.exp.stimID));
         disp(['Added visID to Exp ' num2str(ind)]);
     end
     if all(All(ind).out.exp.visID==0)
         All(ind).out.exp.visID = ones(size(All(ind).out.exp.visID));
         disp(['Corrected VisID to ones'])
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
     
      % create a trial by trial stimSuccess limit
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.stimID);
    numTrials = size(All(ind).out.exp.zdfData,3);
    
    rdata = All(ind).out.exp.rdData;
    bdata = All(ind).out.exp.bdata;
    

    
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
            htg = All(ind).out.exp.holoTargets{h};
            htg(isnan(htg))=[];
            
            vals = rdata(htg,k) - bdata(htg,k);
            stimScore = vals>stimsuccessZ;
            stimSuccessTrial(k)= mean(stimScore) > stimEnsSuccess;
            tempStimScore(k) = mean(stimScore); 
        end
    end
    
    All(ind).out.exp.stimSuccessTrial = stimSuccessTrial;
    percentSuccessStim(ind)=mean(stimSuccessTrial);
    
    clear ensStimScoreExp
    for k=1:numel(us)
        s=us(k);
        ensStimScoreExp(k) = mean(stimSuccessTrial(All(ind).out.exp.stimID==s));
    end
    
    All(ind).out.exp.ensStimScore = ensStimScoreExp;
    ensStimScore{ind}=ensStimScoreExp;
    
     
     
     fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

 end
 

 
 %%Get the number of spikes in each stimulus

clear numSpikesEachStim numCellsEachEns hzEachEns
for ind = 1:numExps
    temp = All(ind).out.exp.stimParams.numPulse;
    numSpikes=[];
    c=0;
    for i=1:numel(temp) %overly complicated way of aligning 0s to be safe if we have 0s that aren't in the begining
        if temp(i)==0
            numSpikes(i)=0;
        else
            c=c+1;
            numSpikes(i) = temp(i)*All(ind).out.exp.stimParams.numCells(c);
        end
    end
    
    
    All(ind).out.anal.numSpikesAddedPerCond = numSpikes;
    numSpikesEachStim{ind} = numSpikes;
    numCellsEachEns{ind} = All(ind).out.exp.stimParams.numCells;
    hzEachEns{ind} = All(ind).out.exp.stimParams.Hz;
    
end
numSpikesEachStim=cell2mat(numSpikesEachStim(:)');
numSpikesEachEns = numSpikesEachStim;
numSpikesEachEns(numSpikesEachStim==0)=[];

numCellsEachEns=cell2mat(numCellsEachEns(:)');
    
hzEachEns = cell2mat(hzEachEns(:)');

ensStimScore=cell2mat(ensStimScore(:)');
ensStimScore(numSpikesEachStim==0)=[];

outVars.ensStimScore=ensStimScore;
outVars.hzEachEns=hzEachEns;
outVars.numCellsEachEns=numCellsEachEns;
outVars.numSpikesEachStim=numSpikesEachStim;
outVars.numSpikesEachEns = numSpikesEachEns;
outVars.percentLowRunTrials=percentLowRunTrials;
