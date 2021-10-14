function [All,outVars] = cleanData(All,opts)
%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=opts.FRDefault;%6;
recWinRange =opts.recWinRange;% [0.5 1.5];% %from vis Start [1.25 2.5];

if isfield(opts,'runValPercent') %back compatible check, moving runVal Percent out of function
    runValPercent = opts.runValPercent;
else
    runValPercent = 0.75;
end
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
    
    if isfield(All(ind).out.exp, 'dataToUse')
        dataToUse = All(ind).out.exp.dataToUse;
    else
        disp(['ind ' num2str(ind) '. no data to use, using zdfData']);
        dataToUse = All(ind).out.exp.zdfData;
    end
    
    All(ind).out.anal.numCells = size(dataToUse,1);
    numCells(ind) = size(dataToUse,1);
    
    if ~isfield(All(ind).out.info,'FR')
        All(ind).out.info.FR=FRDefault;
    end
    
    sz = size(dataToUse);
    
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
    bwinToUse = max(floor([0 All(ind).out.exp.visStart]*All(ind).out.info.FR),[1 1]);
    
    rdata = squeeze(mean(dataToUse(:,winToUse(1):winToUse(2),:),2));
    bdata = squeeze(mean(dataToUse(:,bwinToUse(1):bwinToUse(2),:),2));
    
    All(ind).out.exp.rdData=rdata;
    All(ind).out.exp.bdata=bdata;
    
    All(ind).out.anal.recWinUsed = winToUse;
    All(ind).out.anal.bwinToUse = bwinToUse;
    All(ind).out.anal.recStartTime = visStart;
    All(ind).out.anal.recStartFrame = round(visStart*All(ind).out.info.FR);
    
    
    %runProcessing Section
    runVal = All(ind).out.exp.runVal;
    rnSz = size(runVal);
    runperiod = [1:min(All(ind).out.anal.recWinUsed(2),rnSz(2))];
    
    lowRunVals = mean((runVal(:,runperiod)<runThreshold)');
    
    lowRunTrials = lowRunVals>runValPercent; %percent of frames that need to be below run threshold
    if numel(lowRunTrials)>numel(All(ind).out.exp.stimID)
        lowRunTrials = lowRunTrials(1:numel(All(ind).out.exp.stimID));
    end
    All(ind).out.exp.lowRunTrials = lowRunTrials;
    
    percentLowRunTrials(ind) = mean(lowRunTrials);
    
    
    %Vis Section
    if isfield(All(ind).out,'vis')
        %         sz2 = size(All(ind).out.vis.zdfData);
        try
            visStart = All(ind).out.vis.visStart;
        catch
            visStart = 0.92; %seemed typical
            All(ind).out.vis.visStart = visStart;
            fprintf('\nError vis.visStart not detected...')
        end
        %         winToUse = min(round((visStart+recWinRange).*All(ind).out.info.FR),[inf sz2(2)]) ;
        %         bwinToUse = max(round([0 visStart]*All(ind).out.info.FR),[1 1]);
        %         %      winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz2(2)]) ;
        %         %      bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);
        %         All(ind).out.anal.visRecWinUsed = winToUse;
        %
        %
        %         rdata = squeeze(mean(All(ind).out.vis.zdfData(:,winToUse(1):winToUse(2),:),2));
        %         bdata = squeeze(mean(All(ind).out.vis.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
        %
        %         All(ind).out.vis.rdata=rdata;
        %         All(ind).out.vis.bdata=bdata;
        
        
        
        
        
        
        if ~isfield(All(ind).out.vis,'lowRunTrials')
            try
                runVal = All(ind).out.vis.runVal;
                rnSz = size(runVal);
                runperiod = [1:min(All(ind).out.anal.recWinUsed(2),rnSz(2))];
                lowRunVals = mean((runVal(:,runperiod)<runThreshold)');
                lowRunTrials = lowRunVals>runValPercent; %percent of frames that need to be below run threshold
                All(ind).out.vis.lowRunTrials = lowRunTrials;
            catch
                disp('Vis Run Trial Problem')
            end
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
    vs = unique(All(ind).out.exp.visID);
    numTrials = size(dataToUse,3);
    
    if isfield(opts, 'stimSuccessByZ') && opts.stimSuccessByZ
        rdata = squeeze(mean(All(ind).out.exp.zdfData(:,winToUse(1):winToUse(2),:),2));
        bdata = squeeze(mean(All(ind).out.exp.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
    else
        rdata = All(ind).out.exp.rdData;
        bdata = All(ind).out.exp.bdata;
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
    
    clear ensStimScoreExp
    for k=1:numel(us)
        s=us(k);
        ensStimScoreExp(k) = mean(stimSuccessTrial(All(ind).out.exp.stimID==s));
    end
    
    All(ind).out.exp.ensStimScore = ensStimScoreExp;
    ensStimScore{ind}=ensStimScoreExp;
    
    
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
    
end

%check if stimParams.roi empty
for ind = 1:numExps
    roiList=[];
    if isfield(All(ind).out.exp.stimParams, 'numCells')
        roiList = [All(ind).out.exp.stimParams.roi{:}];
        roiList(roiList==0) = [];
    end
    
    if isempty(roiList)
        disp(['***Error no rois listed in stimParams.rois ind: ' num2str(ind)]);
    end
end

%check that stimParams.numCells is there
for ind = 1:numExps
    if ~isfield(All(ind).out.exp.stimParams, 'numCells')...
            | numel(All(ind).out.exp.stimParams.Seq)-1 ~= numel(All(ind).out.exp.stimParams.numCells)
        try
            roiList = [All(ind).out.exp.stimParams.roi{:}];
            roiList(roiList==0) = [];
            roiCellCount = cellfun(@(x) numel(x), All(ind).out.exp.rois);
            All(ind).out.exp.stimParams.numCells = roiCellCount(roiList);
            
            disp(['Corrected numCells on Ind: ' num2str(ind)]);
        catch
            disp(['*** ERROR COULD NOT CORRECT NUMCELL IND: ' num2str(ind)]);
        end
    end
end

%check that Hz exists
for ind = 1:numExps
    if ~isfield(All(ind).out.exp.stimParams, 'Hz')
        All(ind).out.exp.stimParams.Hz = ones(size(All(ind).out.exp.stimParams.numCells)) * 30;
        disp(['Corrected stimRate to 30Hz (note if incorrect) on ind: ' num2str(ind)])
    end
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
            try
                numSpikes(i) = temp(i)*All(ind).out.exp.stimParams.numCells(c);
            catch
                disp(['Error ind ' num2str(ind) ' stimParams is wrong. assuming all trial types identical (line 261)'])
                if i>1
                    numSpikes(i) = numSpikes(i-1);
                else
                    numSpikes(i) = 0
                end
            end
        end
    end
    
    
    All(ind).out.anal.numSpikesAddedPerCond = numSpikes;
    numSpikesEachStim{ind} = numSpikes;
    numCellsEachEns{ind} = All(ind).out.exp.stimParams.numCells;
    hzEachEns{ind} = All(ind).out.exp.stimParams.Hz;
    
    if numel( All(ind).out.exp.stimParams.numCells) ~= numel(All(ind).out.exp.stimParams.Hz)
        disp(['Error on the number of items ind: ' num2str(ind)])
    end
    
end
numSpikesEachStim=cell2mat(numSpikesEachStim(:)');
numSpikesEachEns = numSpikesEachStim;
numSpikesEachEns(numSpikesEachStim==0)=[];

numCellsEachEns=cell2mat(numCellsEachEns(:)');
hzEachEns = cell2mat(hzEachEns(:)');
try
    numSpikesEachCell = numSpikesEachEns/numCellsEachEns;
catch
    numSpikesEachCell=0;
end
ensStimScore=cell2mat(ensStimScore(:)');
ensStimScore(numSpikesEachStim==0)=[];

outVars.ensStimScore = ensStimScore;
outVars.hzEachEns = hzEachEns;
outVars.numCellsEachEns = numCellsEachEns;
outVars.numSpikesEachStim = numSpikesEachStim;
outVars.numSpikesEachEns = numSpikesEachEns;
outVars.percentLowRunTrials = percentLowRunTrials;
outVars.numSpikesEachCell = numSpikesEachCell;

%% catch some errors in viscode
disp('Running Additional corrections...')
for ind = 1:numExps
    vs = unique(All(ind).out.exp.visID);
    if ~all(floor(vs)==vs)
        disp(['NonStandard visID... Correcting... ind: ' num2str(ind)]);
        tempVisID = All(ind).out.exp.visID;
        All(ind).out.exp.visIDBackup = tempVisID;
        vs=unique(tempVisID);
        newVisID = zeros(size(tempVisID));
        for i =1:numel(tempVisID)
            newVisID(i) = find(vs==tempVisID(i));
        end
        All(ind).out.exp.visID=newVisID;
    end
end
disp('Done')




