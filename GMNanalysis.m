%% Load Experiments

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

% loadPath = 'U:\ioldenburg\outputdata1'
% loadPath = 'C:\Users\ian\Dropbox\Adesnik\Data\outputdata1'
% loadPath = 'C:\Users\SabatiniLab\Dropbox\Adesnik\Data\outputdata1' %Ian Desktop
% loadPath = 'C:\Users\Will\Local Data\100spikes-results\outfiles-all'
%%
numExps = numel(loadList);
if numExps ~= 0
clear All
if ~iscell(loadList)
    numExps=1;
    temp = loadList;
    clear loadList;
    loadList{1} = temp;
end
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end
else
    disp('Did you press this by accident?')
end

%% Hardcode error fixes
 All(1).out.exp.visIDBackup=All(1).out.exp.visID;

%%and do it
visID = All(1).out.exp.visIDBackup;
visID(visID==2)=6;
All(1).out.exp.visID=visID;

% experiment 11 has all 0 (grey screen)
All(11).out.exp.visID = ones(1,length(All(11).out.exp.visID));
%% experiment 12 somehow has a 30 and 3 in the begining of stimParams
badID = All(12).out.exp.stimID==16;
All(12).out.exp.stimID(badID)=[];
All(12).out.exp.zdfData(:,:,badID)=[];
All(12).out.exp.allData(:,:,badID)=[];
All(12).out.exp.runVal(badID,:)=[];
All(12).out.exp.lowMotionTrials(badID)=[];
All(12).out.exp.visID(badID)=[];
All(12).out.exp.outputsInfo.OutputStims(2)=[];
All(12).out.exp.outputsInfo.OutputNames(2)=[];
All(12).out.exp.outputsInfo.OutputOrder(2)=[];
All(12).out.exp.outputsInfo.OutputPatterns(2)=[];
All(12).out.exp.stimParams.Seq(2)=[];
All(12).out.exp.stimParams.numPulse(2)=[];
All(12).out.exp.stimParams.roi(2)=[];





%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=6;
recWinRange = [0.5 1.5];% %from vis Start [1.25 2.5];


%Stim Success Thresholds
stimsuccessZ = 0.5; %over this number is a succesfull stim
stimEnsSuccess = 0.5; %fraction of ensemble that needs to be succsfull

%run Threshold
runThreshold = 6 ; %trials with runspeed below this will be excluded

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
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);

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
    
    if ind==11
        runVal = All(ind).out.vis.runVal;
        rnSz = size(runVal);
        runperiod = [1:min(All(ind).out.anal.recWinUsed(2),rnSz(2))];
        lowRunVals = mean((runVal(:,runperiod)<runThreshold)');
        lowRunTrials = lowRunVals>0.75; %percent of frames that need to be below run threshold
        All(ind).out.vis.lowRunTrials = lowRunTrials;
    end
    
    
    
    percentLowRunTrials(ind) = mean(lowRunTrials);
    
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

%% Make all dataPlots into matrixes of mean responses
 %%Determine Vis Responsive and Process Correlation
 visAlpha = 0.05;
 
 %oftarget risk params
 thisPlaneTolerance = 11.25;%7.5;%1FWHM%10; %in um;% pixels
 onePlaneTolerance = 22.5;%15;%2FWHM %20;
 



clear popResponse pVisR pVisT
ensIndNumber=[];
for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);

    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;
    
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);

    clear respMat baseMat %Order stims,vis,cells
    for i=1:numel(us)
        s = us(i);

        for k= 1 : numel(vs)
            v = vs(k);
            if v~=0
                respMat(i,v,:) = mean(All(ind).out.exp.rdData(:,...
                    trialsToUse & All(ind).out.exp.stimID ==s &...
                    All(ind).out.exp.visID ==v), 2) ;
                baseMat(i,v,:) = mean(All(ind).out.exp.bdata(:,...
                    trialsToUse & All(ind).out.exp.stimID ==s &...
                    All(ind).out.exp.visID ==v), 2) ;
            end
        end
    end

    All(ind).out.anal.respMat = respMat;
    All(ind).out.anal.baseMat = baseMat;


    %%offtargetRisk
    stimCoM = All(ind).out.exp.stimCoM;
    numCells = size(All(ind).out.exp.zdfData,1);
    allCoM = All(ind).out.exp.allCoM;
    stimDepth = All(ind).out.exp.stimDepth;
    allDepth = All(ind).out.exp.allDepth;
    muPerPx = 800/512;

    allLoc = [allCoM*muPerPx (allDepth-1)*30];
    stimLoc = [stimCoM*muPerPx (stimDepth-1)*30];

    roisTargets = All(ind).out.exp.rois;
    holoTargets = All(ind).out.exp.holoTargets;

 
    radialDistToStim=zeros([size(stimCoM,1) numCells]);
    axialDistToStim = zeros([size(stimCoM,1) numCells]);
    StimDistance = zeros([size(stimCoM,1) numCells]);
    for i=1:size(stimCoM,1);
        for k=1:numCells;
            D = sqrt(sum((stimCoM(i,:)-allCoM(k,:)).^2));
            radialDistToStim(i,k)=D;
            z = stimDepth(i)-allDepth(k);
            axialDistToStim(i,k) = z;
            StimDistance(i,k) = sqrt(sum((stimLoc(i,:)-allLoc(k,:)).^2));

        end
    end

    offTargetRisk = zeros([numel(roisTargets) numCells]);
    for i=1:numel(roisTargets)
        Tg = roisTargets{i};
        try
        TgCells = holoTargets{i};
        catch;end;
        
        if numel(Tg) == 1
            temp = radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0;
            temp2 = radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1;
        else
            temp = any(radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0);
            temp2 = any(radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1);
        end
        offTargetRisk(i,:) = temp | temp2;
    end
    All(ind).out.anal.offTargetRisk = offTargetRisk;


    %%ROIinArtifact
    try
        yoffset = -All(ind).out.info.offsets(2);
    catch
        yoffset = 0 ;
    end

    ArtifactSizeLeft = 100;
    ArtifactSizeRight = 100;
    ROIinArtifact = allCoM(:,2)<ArtifactSizeLeft-yoffset | allCoM(:,2)>511-(ArtifactSizeRight+yoffset);
    All(ind).out.anal.ROIinArtifact = ROIinArtifact;
%     pVisR = All(ind).out.anal.pVisR;
%     pVisT = All(ind).out.anal.pVisT;
    


    %ID tuned Cells, should comparing no contrast to with contrast
    pVisR=[];%pVisT=[];
    for i=1:All(ind).out.anal.numCells
        trialsToUse = All(ind).out.exp.visID~=0 &...
            All(ind).out.exp.lowMotionTrials &...
            All(ind).out.exp.lowRunTrials &...
            All(ind).out.exp.stimID==min(All(ind).out.exp.stimID);
 %(All(ind).out.exp.visID==1 | All(ind).out.exp.visID==max(All(ind).out.exp.visID) ) & All(ind).out.exp.lowMotionTrials;



        pVisR(i) = anova1(All(ind).out.exp.rdData(i,trialsToUse),All(ind).out.exp.visID(trialsToUse),'off');
%          pVisR(i) = ranksum(All(ind).out.exp.rdData(i,trialsToUse & All(ind).out.exp.visID==1),...
%              All(ind).out.exp.rdData(i,trialsToUse & All(ind).out.exp.visID== max(All(ind).out.exp.visID)) );
         
%         trialsToUse = All(ind).out.vis.visID~=0 & All(ind).out.vis.visID~=1 & All(ind).out.vis.lowMotionTrials;
%         pVisT(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
    end
    % All(ind).out.anal.pVisR = pVisR;
    
     % hack for expt 11
    % use the vis stim expt
    if ind==11
        pVisR=[];
        for i=1:All(ind).out.anal.numCells     
            trialsToUse = All(ind).out.vis.visID~=0 &...
                All(ind).out.vis.lowMotionTrials &...
                All(ind).out.vis.lowRunTrials;
            pVisR(i) = anova1(All(ind).out.vis.rdData(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
        end
    end
    All(ind).out.anal.pVisR = pVisR;
     
    
    All(ind).out.anal.visPercent = sum(pVisR<visAlpha) / numel(pVisR);
    visPercent(ind) =  All(ind).out.anal.visPercent;

    %%Get Pop Responses
    %         v=1; %best bet for no vis stim.
    vs(vs==0)=[];
    clear popResp popRespDist popRespDistNumCells popRespDistSubtracted  popRespDistVisNumCells popRespDistSubVis popRespDistVis
    clear minDistbyHolo geoDistbyHolo meanDistbyHolo harmDistbyHolo
    
    
    distBins = [0:25:1000];
    numDist =numel(distBins)-1;
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    numVis = numel(vs);

    
    popResp = nan([numStims max(vs)]);
    popRespDist = nan([numStims max(vs) numDist]);
    popRespDistNumCells = nan([numStims max(vs) numDist]);
    popRespDistSubtracted = nan([numStims max(vs) numDist]);
    popRespDistSubVis = nan([numStims max(vs) numDist]);
    popRespDistVis = nan([numStims max(vs) numDist]);
    popRespDistVisNumCells = nan([numStims max(vs) numDist]);
    for k = 1:numVis
        v=vs(k); 
        for i= 1:numStims
            holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
            if i==1;
                cellsToUse = ~ROIinArtifact';
            else
                cellsToUse = ~ROIinArtifact'  & ~offTargetRisk(holo,:);
            end
            popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if i~=1
                Tg=All(ind).out.exp.rois{holo};
                dists = StimDistance(Tg,:);
                
                minDist = min(dists,[],1);
                geoDist = geomean(dists,1); 
                meanDist = mean(dists,1);
                harmDist = harmmean(dists,1);
                
                minDistbyHolo(i,:) = minDist;
                geoDistbyHolo(i,:) = geoDist;
                meanDistbyHolo(i,:) = meanDist;
                harmDistbyHolo(i,:) = harmDist;
                
                distToUse = minDist; % CHANGE THIS (when you want to change whats being analyzed)

%                 cellsToUse = ~ROIinArtifact' &...
%                         ~offTargetRisk(holo,:) &...
%                         minDist > 75; 
%                  popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                
%                 distBins = [0:25:1000]; %moved up
                for d = 1:numel(distBins)-1
                    cellsToUse = ~ROIinArtifact' &...
                        ~offTargetRisk(holo,:) &...
                        distToUse > distBins(d) &...
                        distToUse <= distBins(d+1) ;
                    popRespDist(i,v,d) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                    popRespDistNumCells(i,v,d) = sum(cellsToUse);
                    noHoloPopResponse = mean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse)));
                    popRespDistSubtracted(i,v,d) = popRespDist(i,v,d) - noHoloPopResponse;
                    
                    cellsToUse = cellsToUse & pVisR<visAlpha;
                    tempResp = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                    noHoloPopResponse = mean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse)));
                    popRespDistSubVis(i,v,d) = tempResp - noHoloPopResponse;
                    popRespDistVis(i,v,d) = tempResp;
                    popRespDistVisNumCells(i,v,d) = sum(cellsToUse);

                  
                end

            end
        
        
        end
    end
    
       All(ind).out.anal.minDistbyHolo = minDistbyHolo;
        All(ind).out.anal.geoDistbyHolo = geoDistbyHolo;
        All(ind).out.anal.meanDistbyHolo = meanDistbyHolo;
        
        
    VisCondToUse = 1; %1 is no vis
    if VisCondToUse > size(popResp,2) 
        disp(['VisCond Not available ind: ' num2str(ind)])
        popResponse{ind} = single(nan(size(popResp(:,1))));
        popResponseDist{ind} = single(nan(size(squeeze(popRespDist(:,1,:)))));
        popResponseNumCells{ind} = double(nan(size(squeeze(popRespDistNumCells(:,1,:)))));
    else
        popResponse{ind} = popResp(:,VisCondToUse);
        popResponseDist{ind} = squeeze(popRespDist(:,VisCondToUse,:));
        popResponseDistVis{ind} = squeeze(popRespDistVis(:,VisCondToUse,:));

        popResponseNumCells{ind} = squeeze(popRespDistNumCells(:,VisCondToUse,:));
    end
    popResponseAll{ind} = popResp; %pop Response by Holo
    popResponseAllDist{ind} = popRespDist; %pop Response by Holo and Distance
    popResponseAllDistSub{ind} = popRespDistSubtracted; %pop Response by Holo and Distance with no holostim subtracted aka: holo evoked response
    popResponseAllNumCells{ind} = popRespDistNumCells; %num cells by holo and distance
    
    popResponseAllDistSubVis{ind} = popRespDistSubVis; %pop response by holo and distance with no holo subtracted only from Vis Cells
    popResponseAllDistVis{ind} = popRespDistVis; %Pop response by HOlo and Distance only from Vis Cells
    popResponseAllDistSubVisNC{ind} = popRespDistVisNumCells; %num cells visR by holo and distance
    
    ensIndNumber = [ensIndNumber ones(size(popResp(:,1)'))*ind];
    

    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

popResponse = cell2mat(popResponse(:));
popResponseEns=popResponse;
popResponseEns(numSpikesEachStim==0)=[];

ensIndNumber(numSpikesEachStim==0)=[];

noStimPopResp = popResponse(numSpikesEachStim==0);

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<0.05)); %remove low vis responsive experiments

%% main Ensembles to Use section
% ensemblesToUse = numSpikesEachEns > 75 & numSpikesEachEns <125 & highVisPercentInd & ensIndNumber~=15 & ensIndNumber~=16; %& numCellsEachEns>10 ;

numTrialsPerEns =[];
for ind=1:numExps
    us=unique(All(ind).out.exp.stimID);
    
    for i=1:numel(us)
         trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(i) ;
    
    numTrialsPerEns(end+1)=sum(trialsToUse);
    end
    
end
numTrialsPerEns(numSpikesEachStim==0)=[];


highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<0.05)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));

excludeInds = ismember(ensIndNumber,[2]); %Its possible that the visStimIDs got messed up




ensemblesToUse = numSpikesEachEns > 75 &...
    numSpikesEachEns <110 &...
    ...highVisPercentInd &...
    lowRunInds &...
    ensStimScore > 0.25 &... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    ~excludeInds &...
    numTrialsPerEns > 10;%&...  
     %& numCellsEachEns>10 ;

indsSub = ensIndNumber(ensemblesToUse);
IndsUsed = unique(ensIndNumber(ensemblesToUse));

sum(ensemblesToUse)
%% Create time series plot
minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));

clear allMeanTS allStdTS allnumTS allMeanTSVis allStdTSVis allnumTSVis
for ind=1:numExps
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    strtFrame = All(ind).out.anal.recStartFrame;
    newStart = strtFrame-minStrtFrame+1;
    
    pVisR = All(ind).out.anal.pVisR;
    
    clear mRespTS sRespTS nResp
    for i = 1:numel(us)
        s = us(i);
        h = All(ind).out.exp.stimParams.roi{i};
        
        if h>0
            tg = All(ind).out.exp.holoTargets{h};
            tg(isnan(tg))=[];
        else
            tg=[];
        end
        cellList = 1:numel(ROIinArtifact);
        
        if i==1
            cellsToUse = ~ROIinArtifact';% & pVisR<0.05 ;
        else
            cellsToUse = ~ROIinArtifact' &...
                ~offTargetRisk(h,:) &...
                ~ismember(cellList,tg);% & pVisR<0.05;
        end
        
        for k=1:numel(vs)
            v=vs(k);
            
            dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end,trialsToUse &...
                All(ind).out.exp.stimID==s &...
                All(ind).out.exp.visID==v );
            
            mDat = mean(dat,3);
            mmDat = mean(mDat,1); %pop Average
            sdDat = std(mDat);
            nDat = size(mDat,1);
            
            
            mRespTS(i,:,k) = mmDat; % mean response time series
            sRespTS(i,:,k) = sdDat; % std response time series (by cell);
            nResp(i,k) = nDat;
        end
        
    end
    All(ind).out.anal.mRespTS= mRespTS;
    All(ind).out.anal.sRespTS= sRespTS;
    All(ind).out.anal.nResp = nResp;
    
    
    allMeanTS{ind} = mRespTS(:,:,1);
    allStdTS{ind} = sRespTS(:,:,1);
    allnumTS{ind} = nResp(:,1);
    
    allMeanTSVis{ind} = mRespTS;
    allStdTSVis{ind} = sRespTS;
    allnumTSVis{ind} = nResp;
    
end

baseline=1;

shortestRec = min(cellfun(@(x) size(x,2),allMeanTS));
temp = cellfun(@(x) x(2:end,1:shortestRec),allMeanTS,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
meanTSSquare = cat(1,temp{:});

%no Resp
temp = cellfun(@(x) x(1,1:shortestRec),allMeanTS,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
meanTSSquareNR = cat(1,temp{:});

temp = cellfun(@(x) x(2:end,1:shortestRec),allStdTS,'uniformoutput',0);
stdTSSquare = cat(1,temp{:});

figure(25);clf
subplot(2,2,2)
imagesc(meanTSSquare(ensemblesToUse,:))
% imagesc(meanTSSquare(:,:))
title('stim')

colormap rdbu
caxis([-0.2 0.2])

clear ax
ax(1)=subplot(2,2,4);
fillPlot(meanTSSquare(ensemblesToUse,:),[],'ci');

subplot(2,2,1)
imagesc(meanTSSquareNR(IndsUsed,:))
title('NoStim')

colormap rdbu
caxis([-0.1 0.1])
ax(2) = subplot(2,2,3);
if numel(IndsUsed)>1
fillPlot(meanTSSquareNR(IndsUsed,:),[],'ci');
else
    p= plot(meanTSSquareNR(IndsUsed,:));
    p.Color = 'k';
    p.LineWidth = 2;
end
linkaxes(ax);

%% Plot
f3 = figure(3);
clf(3)


%scatter(meanOSI(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
scatter(1:sum(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')

% p = polyfit(ensOSI(ensemblesToUse),popResponseEns(ensemblesToUse),1);
% f = polyval(p, ensOSI(ensemblesToUse));
% hold on
% plot(ensOSI(ensemblesToUse), f)
% hold off
xlabel('Order of being done')
ylabel('Population Mean Response')
title('OSIs by Ensemble Size')
set(gcf(),'Name','OSIs by Ensemble Size')
cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';

%% group conditions
% not proud of this but it did prevent me from re-writing a bunch of code

numCellsEachEnsBackup = numCellsEachEns;
%numCellsEachEns = numCellsEachEnsBackup;

numCellsEachEns(numCellsEachEns <=5) = 5;
numCellsEachEns(numCellsEachEns > 10) = 20;



%% more simple, take the means, population response by ensemble size
clear avg err ns ens2plt
f6 = figure(6);
clf(f6)
numEns = numel(unique(numCellsEachEns(ensemblesToUse)));
uniqueEns = unique(numCellsEachEns(ensemblesToUse));

x = 1:numEns;
clear data names
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{i} = popResponseEns(ens2plot);
    names{i} = string(uniqueEns(i));
    avg(i) = mean(popResponseEns(ens2plot));
    err(i) = sem(popResponseEns(ens2plot));
    ns(i) = numel(popResponseEns(ens2plot));
end

data{end+1} = noStimPopResp;
names{end+1} = 'No Stim';

cmap=colormap(viridis(numEns));
cmap(end+1,:)=rgb('grey');
p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
% bar(x, avg)
% hold on
% er = errorbar(x, avg, err);
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% hold off
ylabel('Population Response (vis responsive)')
% xticklabels(uniqueEns)
% xticks = 1:6;
title('Mean population response to holo')
xlabel('Ensemble Size')
set(gcf(),'Name','Mean population response to holo')
% ns

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

 r = refline(0);
    r.LineStyle=':';
    r.Color = rgb('grey');

pValEnselbeSize = anovan(popResponseEns(ensemblesToUse),numCellsEachEns(ensemblesToUse)','display','off')
% 
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))

for i=1:size(data,2)
    %     prs = ranksum(data{i},0);
    psr = signrank(data{i});
    
    [h p ] = ttest(data{i},0);
    disp(['Signed Rank: ' num2str(psr,3) '. ttest: ' num2str(p,3)])
end

%% look at just the 10s data for each mouse


% allens2plt = popResponseEns(numCellsEachEns(ensemblesToUse))';

f7 = figure(7);
clf(f7)
k=0;
ens_ids = ensIndNumber(ensemblesToUse);
ens_sizes = numCellsEachEns(ensemblesToUse);
popResponseClip = popResponseEns(ensemblesToUse); %indexing error need to subselect first 

clear sp
for s=unique(ens_sizes)
    clear ens2plt expid exp2plt names
    k=k+1;
    hold on
    sp(k) = subplot(1,numel(unique(numCellsEachEns(ensemblesToUse))),k);
    
    expid = ens_ids(ens_sizes==s);
    ens2plt = popResponseClip(ens_sizes==s)'; %indexing error need to subselect first 

    c=0;
    for i=unique(expid)
        
        c = c+1;
        exp2plt{c} = ens2plt(expid==i);
        names{c}=strrep(All(i).out.info.mouse, '_', '.');
    end

    cmap=colormap(viridis(numel(exp2plt)));
    p=plotSpread(exp2plt,'xNames',names,'showMM',4,'distributionColors',cmap);
    ax=p{3};
    set(findall(gcf(),'type','line'),'markerSize',16)
    p{2}(1).Color = rgb('darkgrey');
    p{2}(2).Color = rgb('darkgrey');
    p{2}(1).LineWidth = 1.5;
    p{2}(2).LineWidth = 1.5;
    uistack(p{2},'bottom')
    xtickangle(45)
    title(['Ensembles of ' num2str(s)])
    
    r = refline(0);
    r.LineStyle=':';
    r.Color = rgb('grey');
    
end

linkaxes(sp(:), 'y')
ax = findobj(sp(1), 'type', 'axes');
set([ax.YLabel], 'string', 'Population Response')
set(gcf(),'Name','Population response to holo by expt and size')
sgtitle('Population response to holo by expt and size')


%% Plot Pop Response by Distance
popDistAll = cell2mat(popResponseDist');
% popDistAll = cell2mat(popResponseDistVis');

popDist = popDistAll(numSpikesEachStim~=0,:);

popNumCells = cell2mat(popResponseNumCells');
popNCells = popNumCells(numSpikesEachStim~=0,:);


ensSizes = unique(numCellsEachEns(ensemblesToUse))   ;


colorList = {rgb('DarkBlue') rgb('steelblue') rgb('gold')};

figure(9);clf
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDist(ensemblesToUse & numCellsEachEns==ensSizes(i),:);
meanDat = nanmean(dat);
stdDat = nanstd(dat);
numpDat = sum(~isnan(dat));
semDat = stdDat./sqrt(numpDat);


hold on
errorbar(distBins(2:end),meanDat,semDat,'linewidth',2,'color',colorList{i})
end
r = refline(0);
r.LineStyle=':';
r.Color = rgb('grey');
r.LineWidth = 2;
xlabel('Distance from a target')
ylabel('Population Response (mean of ensembles'' pop response)')
 xlim([0 350])
legend('Small', 'Medium', 'Big')

% pValEnselbeSize = anovan(popDist(ensemblesToUse,1:20),numCellsEachEns(ensemblesToUse)','display','off')


%% as above but for full and no vis Conditions

figure(11);clf


    
popDatNoVis = cell2mat(cellfun(@(x) squeeze(x(2:end,1,:)), popResponseAllDist,'uniformoutput',0)');
% popDatNoVisNoStim = squeeze(cell2mat(cellfun(@(x) (x(1,1,:)), popResponseAllDist,'uniformoutput',0)'));

popDatMaxVis = cell2mat(cellfun(@(x) squeeze(x(2:end,end,:)), popResponseAllDist,'uniformoutput',0)');
% popDatMaxVisNoStim = squeeze(cell2mat(cellfun(@(x) (x(1,end,:)), popResponseAllDist,'uniformoutput',0)'));
popDatMaxVisSubtracted = cell2mat(cellfun(@(x) squeeze(x(2:end,end,:)), popResponseAllDistSub,'uniformoutput',0)');
popDatMaxVisSubVis = cell2mat(cellfun(@(x) squeeze(x(2:end,end,:)), popResponseAllDistSubVis,'uniformoutput',0)');


divider = 1;
popDatVis2 = cell2mat(cellfun(@(x) squeeze(x(2:end,round(size(x,2)/divider),:)), popResponseAllDist,'uniformoutput',0)');
popDatVis2Subtracted = cell2mat(cellfun(@(x) squeeze(x(2:end,round(size(x,2)/divider),:)), popResponseAllDistSub,'uniformoutput',0)');


ensSizes = unique(numCellsEachEns(ensemblesToUse))   ;
colorList = {rgb('DarkBlue') rgb('steelblue') rgb('gold')};



subplot(1,2,1)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatNoVis(ensemblesToUse & numCellsEachEns==ensSizes(i) ,:);
meanDat = nanmean(dat);
stdDat = nanstd(dat);
numpDat = sum(~isnan(dat));
semDat = stdDat./sqrt(numpDat);
hold on
errorbar(distBins(2:end),meanDat,semDat,'linewidth',2,'color',colorList{i})
end
r = refline(0);
r.LineStyle=':';
r.Color = rgb('grey');
r.LineWidth = 2;
xlabel('Minimal distance from a target')
ylabel('Population Response (mean of ensembles'' pop response)')
xlim([0 350])
legend('Small', 'Medium', 'Big')
title('Holographic induced changes 0 Contrast')

subplot(1,2,2)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatMaxVisSubtracted(ensemblesToUse & numCellsEachEns==ensSizes(i) ,:);
% dat = popDatMaxVis(ensemblesToUse & numCellsEachEns==ensSizes(i) & highVisPercentInd ,:);

meanDat = nanmean(dat);
stdDat = nanstd(dat);
numpDat = sum(~isnan(dat));
semDat = stdDat./sqrt(numpDat);
hold on
errorbar(distBins(2:end),meanDat,semDat,'linewidth',2,'color',colorList{i})
end
r = refline(0);
r.LineStyle=':';
r.Color = rgb('grey');
r.LineWidth = 2;
xlabel('Minimal distance from a target')
ylabel('Population Response (mean of ensembles'' pop response)')
xlim([0 350])
legend('Small', 'Medium', 'Big')
title('Holographic induced changes Max Contrast')


figure(12);clf
contrastsToView = [6 3 2 1.5 1.25 1] ;%I know its weird i just wanted to be able to catch times that we were using different contrasts, will work out to 1:6 if there are 6 contrasts; 1:6;
for c=1:numel(contrastsToView)
ax(c) = subplot(1,numel(contrastsToView),c);
divider = contrastsToView(c);
popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistSub,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDist,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistSubVis,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistVis,'uniformoutput',0)');

for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatToPlot(ensemblesToUse & numCellsEachEns==ensSizes(i) & highVisPercentInd ,:);
meanDat = nanmean(dat);
stdDat = nanstd(dat);
numpDat = sum(~isnan(dat));
semDat = stdDat./sqrt(numpDat);
hold on
errorbar(distBins(2:end),meanDat,semDat,'linewidth',2,'color',colorList{i})
end
r = refline(0);
r.LineStyle=':';
r.Color = rgb('grey');
r.LineWidth = 2;
xlabel('Minimal distance from a target')
ylabel('Population Response (mean of ensembles'' pop response)')
xlim([0 400])
legend('Small', 'Medium', 'Big')
title(['Contrast ' num2str(c) ] );



end

linkaxes([ax(:)])

%% Contrast response functions
distBinCRF = [50:25:350];
visAlphaCRF = 0.05;
counter = 0;
clear CRF CRFNS CRFDist CRFNSDist
for ind = 1:numExps
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    pVisR = All(ind).out.anal.pVisR;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    for i=1:numel(us)
        s = us(i);
        h = All(ind).out.exp.stimParams.roi{i};
        
        if h==0
            continue
        end
        cellsToUse = ~ROIinArtifact' & pVisR<visAlphaCRF &...
            ~any(offTargetRisk);%~offTargetRisk(h,:); %~any(offTargetRisk);% (h,:); 
  
        
        counter=counter+1;
        CRF(counter,:) = nanmean(respMat(i,:,cellsToUse) - baseMat(i,:,cellsToUse),3);
        CRFNS(counter,:) = nanmean(respMat(1,:,cellsToUse) - baseMat(1,:,cellsToUse),3);
        
        for d=1:numel(distBinCRF)-1
            cellsToUseDist = cellsToUse &...
                All(ind).out.anal.minDistbyHolo(h+1,:) <=distBinCRF(d+1) &...
                All(ind).out.anal.minDistbyHolo(h+1,:) >distBinCRF(d) ;
            
            CRFDist(counter,:,d) = nanmean(respMat(i,:,cellsToUseDist) - baseMat(i,:,cellsToUseDist),3);
            CRFNSDist(counter,:,d) = nanmean(respMat(1,:,cellsToUseDist) - baseMat(1,:,cellsToUseDist),3);
        
        end
    end
end

%catch missing data
shouldBeNAN = CRF==0 & CRFNS ==0;
CRF(shouldBeNAN) = nan;
CRFNS(shouldBeNAN) = nan;       
CRFDiff = CRF-CRFNS;

shouldBeNAN = CRFDist==0 & CRFNSDist ==0;
CRFDist(shouldBeNAN) = nan;
CRFNSDist(shouldBeNAN) = nan;       
CRFDiffDist = CRFDist-CRFNSDist;

highVisPercentIndMoreStringent = ~ismember(ensIndNumber,find(visPercent<0.1)); %remove low vis responsive experiments

ensemblesToUseThis = ensemblesToUse;% & highVisPercentIndMoreStringent ;%& ensIndNumber==8;% highVisPercentIndMoreStringent ;

figure(13);clf
subplot(1,2,1)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat1 = CRF(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:);
meanDat = nanmean(dat1);
stdDat = nanstd(dat1);
numpDat = sum(~isnan(dat1));
semDat = stdDat./sqrt(numpDat);

data1{i}=dat1;

hold on
errorbar(meanDat,semDat,'linewidth',2,'color',colorList{i})

% dat2 = CRFNS(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:);

dat2 = CRFNS(ensemblesToUseThis ,:);
indsOrigin{i} = ensIndNumber(ensemblesToUseThis & numCellsEachEns==ensSizes(i));
data2{i}=dat2;
meanDat = nanmean(dat2);
stdDat = nanstd(dat2);
numpDat = sum(~isnan(dat2));
semDat = stdDat./sqrt(numpDat);
e2=errorbar(meanDat,semDat,'linewidth',2,'color',rgb('grey'),'linestyle',':');
% e2.Color = colorList{i};

end

ylabel('Z-score dF')
xticks(1:6);
xticklabels({'0%' '1%' '4%' '10%' '40%' '100%'})
xlabel('Contrast')

subplot(1,2,2)

for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat1 = CRFDiff(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:);
meanDat = nanmean(dat1);
stdDat = nanstd(dat1);
numpDat = sum(~isnan(dat1));
semDat = stdDat./sqrt(numpDat);

data1{i}=dat1;

hold on
errorbar(meanDat,semDat,'linewidth',2,'color',colorList{i})
end
r = refline(0);
r.LineWidth = 2;
r.LineStyle =':';
r.Color = rgb('grey');

ylabel('\Delta Z-Scored dF')
xlabel('Contrast')
xticks(1:6);
xticklabels({'0%' '1%' '4%' '10%' '40%' '100%'})

clear ax
figure(14);clf
for c=1:numel(distBinCRF)-1
    subplot(1,numel(distBinCRF)-1,c)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
% dat1 = CRFDiffDist(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:,c);
dat1 = CRFDist(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:,c);

meanDat = nanmean(dat1);
stdDat = nanstd(dat1);
numpDat = sum(~isnan(dat1));
semDat = stdDat./sqrt(numpDat);

data1{i}=dat1;

hold on
errorbar(meanDat,semDat,'linewidth',2,'color',colorList{i})

% dat2 = CRFNSDist(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:,c);

dat2 = CRFNSDist(ensemblesToUseThis ,:,c);
indsOrigin{i} = ensIndNumber(ensemblesToUseThis & numCellsEachEns==ensSizes(i));
data2{i}=dat2;
meanDat = nanmean(dat2);
stdDat = nanstd(dat2);
numpDat = sum(~isnan(dat2));
semDat = stdDat./sqrt(numpDat);
% errorbar(meanDat,semDat,'linewidth',2,'color',rgb('grey'),'linestyle',':')

end


ylabel('Z-score dF')
xticks(1:6);
xticklabels({'0%' '1%' '4%' '10%' '40%' '100%'})
xlabel('Contrast')

title([num2str(distBinCRF(c)) ' to ' num2str(distBinCRF(c+1)) ' \mum from Target'])
end
linkaxes

%% Calculate L1 and L2 
contrastsToView = [6 3 2 1.5 1.25 1] ;%I know its weird i just wanted to be able to catch times that we were using different contrasts, will work out to 1:6 if there are 6 contrasts; 1:6;

clear EnsL1 EnsL2
c=0;
for ind =1:numExps
    
  for h= 1:numel(All(ind).out.exp.stimParams.Seq)-1
            holo = All(ind).out.exp.stimParams.roi{h+1}; % only cycle through holos
    
            divider =inf; %1 is max vis. inf is 0
            maxV = max(All(ind).out.exp.visID);
            v = max(round(maxV/divider),1);
            
            trialsToUse = All(ind).out.exp.lowMotionTrials & All(ind).out.exp.visID==v;
%             cellsToUse =  ~All(ind).out.anal.ROIinArtifact' & All(ind).out.anal.offTargetRisk(holo,:);
                        cellsToUse =  ~All(ind).out.anal.ROIinArtifact' & ~any(All(ind).out.anal.offTargetRisk(:,:));

                        
            us = unique(All(ind).out.exp.stimID); 
            
            testData = All(ind).out.exp.rdData(cellsToUse,trialsToUse & All(ind).out.exp.stimID == us(h+1));
            ExpectedData = All(ind).out.exp.rdData(cellsToUse,trialsToUse & All(ind).out.exp.stimID == us(1));
            [L1 L2 L3] =  calcL1L2(testData,ExpectedData,2);
            L1= L1/size(testData,1);
            L2 = L2/sqrt(size(testData,1));
            L3 = L3/ ((size(testData,1))^(1/3)) ;
            
            
            c=c+1;
            EnsL1(c) = L1;
            EnsL2(c) = L2;      
  end
    
end


clear avg err ns ens2plt
f6 = figure(10);
clf(f6)
numEns = numel(unique(numCellsEachEns(ensemblesToUse)));
uniqueEns = unique(numCellsEachEns(ensemblesToUse));

subplot(1,3,1)
x = 1:numEns;
clear data names
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{i} = EnsL1(ens2plot);
    names{i} = string(uniqueEns(i));
    avg(i) = mean(popResponseEns(ens2plot));
    err(i) = sem(popResponseEns(ens2plot));
    ns(i) = numel(popResponseEns(ens2plot));
end

% data{end+1} = noStimPopResp;
% names{end+1} = 'No Stim';

cmap=colormap(viridis(numEns));
% cmap(end+1,:)=rgb('grey');
p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
ylabel('L1')
% xticklabels(uniqueEns)
% xticks = 1:6;
% title('Mean population response to holo')
xlabel('Ensemble Size')
set(gcf(),'Name','L1 of Nonstimulated Cells to holo')

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

subplot(1,3,2)
x = 1:numEns;
clear data names
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{i} = EnsL2(ens2plot);
    names{i} = string(uniqueEns(i));
    avg(i) = mean(popResponseEns(ens2plot));
    err(i) = sem(popResponseEns(ens2plot));
    ns(i) = numel(popResponseEns(ens2plot));
end

% data{end+1} = noStimPopResp;
% names{end+1} = 'No Stim';

cmap=colormap(viridis(numEns));
% cmap(end+1,:)=rgb('grey');
p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
ylabel('L2')
% xticklabels(uniqueEns)
% xticks = 1:6;
% title('Mean population response to holo')
xlabel('Ensemble Size')
set(gcf(),'Name','L2 of Nonstimulated Cells to holo')

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

subplot(1,3,3)
x = 1:numEns;
clear data names group
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{i} = EnsL1(ens2plot)./EnsL2(ens2plot);
    group{i} = ones([1 numel(ens2plot)])*uniqueEns(i);
    names{i} = string(uniqueEns(i));
    avg(i) = mean(popResponseEns(ens2plot));
    err(i) = sem(popResponseEns(ens2plot));
    ns(i) = numel(popResponseEns(ens2plot));
end

% data{end+1} = noStimPopResp;
% names{end+1} = 'No Stim';

cmap=colormap(viridis(numEns));
% cmap(end+1,:)=rgb('grey');
p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
ylabel('L1/L2')
ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

pValEnselbeSize = anova1([data{:}]',[group{:}]','off')
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))

% for i=1:size(data,2)
%     %     prs = ranksum(data{i},0);
%     psr = signrank(data{i});
%     
%     [h p ] = ttest(data{i},0);
%     disp(['Signed Rank: ' num2str(psr,3) '. ttest: ' num2str(p,3)])
% end
%% L1 L2 by dist
All(ind).out.anal.minDistbyHolo;
distBins = [0:25:500];

calcMethod=2;

clear EnsL1 EnsL2
c=0;
for ind =1:numExps
    
  for h= 1:numel(All(ind).out.exp.stimParams.Seq)-1
            holo = All(ind).out.exp.stimParams.roi{h+1}; % only cycle through holos
%             divider = 1; %6 is no vis, 1 is max vis
%             x = max(All(ind).out.exp.visID);
%             v = max(round(size(x,2)/divider),1);
            numTrials = numel(All(ind).out.exp.stimID);
            trialList = 1:numTrials;

            trialsToUse = All(ind).out.exp.lowMotionTrials &...
                All(ind).out.exp.lowRunTrials &...
                All(ind).out.exp.stimSuccessTrial &...
                All(ind).out.exp.visID==1 ;%&...
%                 isEven(trialList);
            cellsToUse =  ~All(ind).out.anal.ROIinArtifact' & ~All(ind).out.anal.offTargetRisk(holo,:);
            
            us = unique(All(ind).out.exp.stimID); 
            
            c=c+1;
            
            for d = 1:numel(distBins)-1
                D = distBins(d+1);
                cellsToUseDist = cellsToUse &...
                    All(ind).out.anal.minDistbyHolo(h+1,:) <=distBins(d+1) &...
                    All(ind).out.anal.minDistbyHolo(h+1,:) >distBins(d) ;
                
                thisTrialsToUse = trialsToUse & All(ind).out.exp.stimID == us(h+1);
%                 if sum(thisTrialsToUse)<3
%                     thisTrialsToUse = logical(zeros(size(thisTrialsToUse)));
%                 end
                tttu(c)=sum(thisTrialsToUse); %it stats for This Trials To Use. its not totaly arbitrary, also i'm sorry. 
                testData = All(ind).out.exp.rdData(cellsToUseDist,thisTrialsToUse);
                ExpectedData = All(ind).out.exp.rdData(cellsToUseDist,trialsToUse & All(ind).out.exp.stimID == us(1));
                [L1 L2] =  calcL1L2(testData,ExpectedData,calcMethod);
                L1= L1/size(testData,1);
                L2 = L2/sqrt(size(testData,1));
                
                EnsL1(c,d) = L1;
                EnsL2(c,d) = L2;
            
                [L1 L2] =  calcL1L2(ExpectedData,ExpectedData,calcMethod);
                L1= L1/size(ExpectedData,1);
                L2 = L2/sqrt(size(ExpectedData,1));
                
                   
                ctrL1(c,d) = L1;
                ctrL2(c,d) = L2;
            end
  end
    
end



figure(111);clf
subplot(1,3,1)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse & tttu>2);
    data = EnsL1(ens2plot,:);
    
    
    e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data)./sqrt(sum(~isnan(data))));
%         e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data));

    e.Color = colorList{i};
    e.LineWidth = 2;
%     names{i} = string(uniqueEns(i));
%     avg(i) = mean(popResponseEns(ens2plot));
%     err(i) = sem(popResponseEns(ens2plot));
%     ns(i) = numel(popResponseEns(ens2plot));

end
ens2plot = find(ensemblesToUse & tttu>2);
ctrData = ctrL1(ens2plot,:);
e2 = errorbar(distBins(1:end-1),nanmean(ctrData,1),nanstd(ctrData)./sqrt(sum(~isnan(ctrData))));
e2.Color = rgb('grey');
e2.LineWidth =2;
e2.LineStyle =':';

% xlim([0 250])
legend('small', 'medium', 'large')
title('L1')
xlabel('Distance')

subplot(1,3,2)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse & tttu>2);
    data = EnsL2(ens2plot,:);
    
    e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data)./sqrt(sum(~isnan(data))));
%         e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data));

    e.Color = colorList{i};
    e.LineWidth = 2;
%     names{i} = string(uniqueEns(i));
%     avg(i) = mean(popResponseEns(ens2plot));
%     err(i) = sem(popResponseEns(ens2plot));
%     ns(i) = numel(popResponseEns(ens2plot));
end
ens2plot = find(ensemblesToUse & tttu>2);
ctrData = ctrL2(ens2plot,:);
e2 = errorbar(distBins(1:end-1),nanmean(ctrData,1),nanstd(ctrData)./sqrt(sum(~isnan(ctrData))));
e2.Color = rgb('grey');
e2.LineWidth =2;
e2.LineStyle =':';

% xlim([0 250])
legend('small', 'medium', 'large')
title('L2')
xlabel('Distance')

subplot(1,3,3)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse & tttu>2);
    data = EnsL2(ens2plot,:)./EnsL1(ens2plot,:);
    
    
    e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data)./sqrt(sum(~isnan(data))));
%         e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data));

    e.Color = colorList{i};
    e.LineWidth = 2;
%     names{i} = string(uniqueEns(i));
%     avg(i) = mean(popResponseEns(ens2plot));
%     err(i) = sem(popResponseEns(ens2plot));
%     ns(i) = numel(popResponseEns(ens2plot));
end

ens2plot = find(ensemblesToUse & tttu>2);
ctrData = ctrL2(ens2plot,:)./ctrL1(ens2plot,:);
e2 = errorbar(distBins(1:end-1),nanmean(ctrData,1),nanstd(ctrData)./sqrt(sum(~isnan(ctrData))));
e2.Color = rgb('grey');
e2.LineWidth =2;
e2.LineStyle =':';
% xlim([0 250])
legend('small', 'medium', 'large')
title('L2/L1 Ratio')
xlabel('Distance')


%% Correlation Analyisis Determine Correlation Coefficients

for ind = 1:numExps
      pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    us = unique(All(ind).out.exp.stimID);

    %Spont Corr - correlation coefficient on time series from no stim
    %period
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) &...
        All(ind).out.exp.visID == 1;
    unrollData = All(ind).out.exp.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    [SpontCorr SpCoP] = corr(unrollData');
    
    %AllCorr - the correlation coef on all time series
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) ;
    unrollData = All(ind).out.exp.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    [AllCorr AlCoP] = corr(unrollData');
    
    %All corr mean - correlation coef of response (not time series)
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) ;
    unrollData = All(ind).out.exp.rdData(:,trialsToUse);
    sz = size(unrollData);
    
    [AllMCorr AmCoP] = corr(unrollData');
    
    %noise corr - correlation coef of residual trial response (not time
    %series) i.e. trial response - mean response for that condition
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) ;
    vs = unique(All(ind).out.exp.visID);
    unrollData = [];
    meanResps = [];
    for k = 1:numel(vs)
        v = vs(k);
        trialsToUseThis = trialsToUse & All(ind).out.exp.visID==v;
        
        dataPart = All(ind).out.exp.rdData(:,trialsToUseThis);
        mData = mean(dataPart');
        meanResps(k,:) =  mData;
        dataPart = dataPart-mData';
        
        unrollData = cat(2,unrollData,dataPart);
    end
        
    [NoiseCorr NoCoP] = corr(unrollData');
    [SignalCorr SiCoP] = corr(meanResps);

    
    All(ind).out.anal.SpontCorr = SpontCorr;
    All(ind).out.anal.SpCoP = SpCoP;

    All(ind).out.anal.AllCorr = AllCorr;
    All(ind).out.anal.AlCoP = AlCoP;
    
    All(ind).out.anal.AllMCorr = AllMCorr;
    All(ind).out.anal.AmCoP =AmCoP;
    
    All(ind).out.anal.SignalCorr = SignalCorr;
    All(ind).out.anal.SiCoP =SiCoP;

    All(ind).out.anal.NoiseCorr = NoiseCorr;
    All(ind).out.anal.NoCoP = NoCoP;

    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%%Determine the Ensemble CoCorrelation

ensSpCo=[];ensAlCo=[];ensAmCo=[];ensSiCo=[];ensNoCo=[];
for ind = 1:numExps
    clear ensembleSoCo ensembleAlCo ensembleAmCo ensembleSiCo ensembleNoCo
    for i =1:numel(All(ind).out.exp.holoTargets)
        ht = All(ind).out.exp.holoTargets{i};
        ht(isnan(ht))=[];
        corrToUse = All(ind).out.anal.SpontCorr;
        corMat = corrToUse(ht,ht); 
        corMat(logical(eye(numel(ht))))=nan;
        ensembleSpCo(i) = nanmean(corMat(:));
        
        corrToUse = All(ind).out.anal.AllCorr;
        corMat = corrToUse(ht,ht); 
        corMat(logical(eye(numel(ht))))=nan;
        ensembleAlCo(i) =nanmean(corMat(:));
        
        corrToUse = All(ind).out.anal.AllMCorr;
        corMat = corrToUse(ht,ht); 
        corMat(logical(eye(numel(ht))))=nan;
        ensembleAmCo(i) = nanmean(corMat(:));
        
        corrToUse = All(ind).out.anal.SignalCorr;
        corMat = corrToUse(ht,ht); 
        corMat(logical(eye(numel(ht))))=nan;
        ensembleSiCo(i) =nanmean(corMat(:));
        
        corrToUse = All(ind).out.anal.NoiseCorr;
        corMat = corrToUse(ht,ht); 
        corMat(logical(eye(numel(ht))))=nan;
        ensembleNoCo(i) =nanmean(corMat(:));
    end
    All(ind).out.anal.ensembleSpCo = ensembleSpCo;
    All(ind).out.anal.ensembleAlCo = ensembleAlCo;
    All(ind).out.anal.ensembleAmCo = ensembleAmCo;
    All(ind).out.anal.ensembleSiCo = ensembleSiCo;
    All(ind).out.anal.ensembleNoCo = ensembleNoCo;
    
    ensSpCo = cat(2,ensSpCo,ensembleSpCo);
    ensAlCo = cat(2,ensAlCo,ensembleAlCo);
    ensAmCo = cat(2,ensAmCo,ensembleAmCo);
    ensSiCo = cat(2,ensSiCo,ensembleSiCo);
    ensNoCo = cat(2,ensNoCo,ensembleNoCo);
end
figure(14);clf
dat = {ensSpCo(ensemblesToUse), ensAlCo(ensemblesToUse), ensAmCo(ensemblesToUse),  ensSiCo(ensemblesToUse), ensNoCo(ensemblesToUse)}; 
names = {'Spont' 'All' 'All (v2)' 'Signal' 'Noise'};
fancyPlotSpread(dat,names);
title('Ensemble Mean Correlations by type')
ylabel('Correlation (Rho)')

%%plot Pop Response by Correlation
f3 = figure(3);
clf(3)

for i=1:5
    subplot(5,1,i)
    dataToUse = dat{i};
scatter(dataToUse,popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
% scatter(1:sum(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')

title([names{i} ' Correlation'])

xlabel(['Correlation of Ensemble'])
ylabel('Population Mean Response')
% title('OSIs by Ensemble Size')
set(gcf(),'Name','OSIs by Ensemble Size')
cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';
end


%% Plot 
clear popResponseCorr
for ind = 1:numExps

    corrToUse  = All(ind).out.anal.AllCorr;

    
    vs =  unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    clear popRespCorr minDistbyHolo cellsToUse
 for v = 1:numel(vs)
        for i= 1:numel(All(ind).out.exp.stimParams.Seq)
            holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
            if i==1
                cellsToUse = ~ROIinArtifact';
            else
                cellsToUse = ~ROIinArtifact'  & ~offTargetRisk(holo,:);
            end
%             popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if i~=1
                Tg=All(ind).out.exp.holoTargets{holo};
                Tg(isnan(Tg))=[];
                
                distCorr = corrToUse(Tg,:);
                minDist = mean(distCorr);
                
                if numel(Tg)==0
                    minDistbyHolo(i,:) = ones([1 size(minDist,2)])*1000;
                else
                    minDistbyHolo(i,:) = minDist;
                end
                distBins = linspace(-0.5,0.5,40);
                for d = 1:numel(distBins)-1
                    cellsToUse = ~ROIinArtifact' &...
                        ~offTargetRisk(holo,:) &...
                        minDist > distBins(d) &...
                        minDist <= distBins(d+1) ;
                    popRespCorr(i,v,d) = nanmean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                    
                    noHoloEquivalent = nanmean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse)));
                    popRespCorrSub(i,v,d) =  popRespCorr(i,v,d) - noHoloEquivalent;
                end
            end
        end
 end
 popRespCorr(1,:,:)=[];
 popRespCorrSub(1,:,:)=[];
 
 popResponseCorr{ind} = popRespCorr;
 popResponseCorrSub{ind} = popRespCorrSub;
end


% 
% popResponseCorr = cell2mat(popResponseCorr(:));
% popResponseCorr(numSpikesEachStim==0)=[];

popRespCorr;


temp = cellfun(@(x) squeeze(x(:,1,:)),popResponseCorr,'uniformoutput',0) ;
% temp = cellfun(@(x) squeeze(x(:,end,:)),popResponseCorr,'uniformoutput',0) ;

% temp = cellfun(@(x) squeeze(x(:,end,:)),popResponseCorrSub,'uniformoutput',0) ;
% temp = cellfun(@(x) squeeze(x(:,1,:)),popResponseCorrSub,'uniformoutput',0) ;

% temp = cellfun(@(x) squeeze(x(:,round(end/2),:)),popResponseCorr,'uniformoutput',0) ;

EnsCorR = cat(1,temp{:});

figure(15);clf
% subplot(1,3,1)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse );
    data = EnsCorR(ens2plot,:);
    
    
    e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data)./sqrt(sum(~isnan(data))));
    
%         e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data));

    e.Color = colorList{i};
    e.LineWidth = 2;
    hack = num2cell(distBins(1:end-1));
    hack = cellfun(@(x) num2str(x),hack,'uniformoutput',0);
%     p = plotSpread(data, 'xNames', hack, 'showMM', 4);
% fancyPlotSpread(data,hack)
%     names{i} = string(uniqueEns(i));
%     avg(i) = mean(popResponseEns(ens2plot));
%     err(i) = sem(popResponseEns(ens2plot));
%     ns(i) = numel(popResponseEns(ens2plot));
end
xlim([-0.4 0.4])
legend('small', 'medium', 'large')
r = refline(0);
r.LineStyle = ':';
r.Color = rgb('grey');
r.LineWidth = 2;

ylabel('Pop Resp to HoloStim')
xlabel('Responder to Ensemble Correlation')

%%

%% Within Ensemble Size Correlation by ensemble OSI
figure(16);clf
hold on

corrToUse = ensAlCo;%meanOSI';

lowCor = 0.025;
highCor = 0.07;

EnsSizeToUse = 10;

temp = cellfun(@(x) reshape(x(:,1,:),size(x(:,1,:),1),size(x(:,1,:),3)),popResponseCorr,'uniformoutput',0) ;
% temp = cellfun(@(x) reshape(x(:,end,:),size(x(:,end,:),1),size(x(:,end,:),3)),popResponseCorr,'uniformoutput',0) ;

EnsCorR = cat(1,temp{:});

CorRange = [0 lowCor highCor 1];
numOfEnsUsed=[];
for i = 1:numel(CorRange)-1
    ens2plot =  ensemblesToUse &...
        numCellsEachEns==EnsSizeToUse &...
        corrToUse >= CorRange(i) &  corrToUse < CorRange(i+1)  ; %numCellsEachEns==EnsSizeToUse &...
    
    data = EnsCorR(ens2plot,:);
    numOfEnsUsed(i) = size(data,1);
    try
        e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data)./sqrt(sum(~isnan(data))));
        e.Color = colorList{i};
        e.LineWidth = 2;
    catch
    end
end
r = refline(0);
r.LineStyle = ':';
r.Color = rgb('grey');
r.LineWidth = 2;
legend(['Low Cor: ' num2str(numOfEnsUsed(1)) ' Ens'],...
    ['Medium Cor: ' num2str(numOfEnsUsed(2)) ' Ens'],...
    ['High Cor: ' num2str(numOfEnsUsed(3)) ' Ens'], 'No Change')


ylabel('Pop Resp to HoloStim')
xlabel('Responder to Ensemble Correlation')
title([num2str(EnsSizeToUse) ' Target Holograms, by OSI'])
