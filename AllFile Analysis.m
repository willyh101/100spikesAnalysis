[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

% loadPath = 'U:\ioldenburg\outputdata1'
% loadPath = 'C:\Users\ian\Dropbox\Adesnik\Data\outputdata1'
% loadPath = 'C:\Users\SabatiniLab\Dropbox\Adesnik\Data\outputdata1' %Ian Desktop

%%
numExps = numel(loadList);
if numExps ~= 0
clear All
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end
else
    disp('Did you press this by accident?')
end

%% Correcting stuff
All(1).out.exp.stimParams.numCells = ones(size(All(1).out.exp.stimParams.numPulse))*20;
All(1).out.exp.stimParams.Hz = ones(size(All(1).out.exp.stimParams.numPulse))*30;

All(2).out.exp.stimParams.numCells = ones(size(All(2).out.exp.stimParams.numPulse))*15;
All(2).out.exp.stimParams.Hz = ones(size(All(2).out.exp.stimParams.numPulse))*30;

All(3).out.exp.stimParams.numCells = ones(size(All(3).out.exp.stimParams.numPulse))*5;
All(3).out.exp.stimParams.Hz = ones(size(All(3).out.exp.stimParams.numPulse))*30;

All(4).out.exp.stimParams.Seq = All(4).out.exp.stimParams.Seq(1:2);

visID = All(23).out.exp.visID;
newVisID = visID;
newVisID(visID==0)=1;
newVisID(visID==1)=2;
newVisID(visID~=0 & visID~=1)=0;
All(23).out.exp.visID = newVisID;


%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=6;
recWinRange = [0.5 1.5];% %from vis Start [1.25 2.5];


%Stim Success Thresholds
stimsuccessZ = 0.5; %over this number is a succesfull stim
stimEnsSuccess = 0.5; %fraction of ensemble that needs to be succsfull

%run Threshold
runThreshold = 6 ; %trials with runspeed below this will be excluded

clear ensStimScore numSpikesEachEns

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
    for i=1:numel(temp); %overly complicated way of aligning 0s to be safe if we have 0s that aren't in the begining
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

    thisPlaneTolerance = 15;10; %in pixels
    onePlaneTolerance = 25;20;

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
     All(ind).out.anal.pVisR = pVisR;
     visAlpha = 0.05;
    
    All(ind).out.anal.visPercent = sum(pVisR<visAlpha) / numel(pVisR);
    visPercent(ind) =  All(ind).out.anal.visPercent;

    %%Get Pop Responses
    %         v=1; %best bet for no vis stim.
    vs(vs==0)=[];
    clear popResp popRespDist popRespDistNumCells popRespDistSubtracted  popRespDistVisNumCells popRespDistSubVis popRespDistVis
    clear minDistbyHolo geoDistbyHolo meanDistbyHolo harmDistbyHolo
    for v = 1:numel(vs)
        for i= 1:numel(All(ind).out.exp.stimParams.Seq)
            holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
            if holo==0;
                cellsToUse = ~ROIinArtifact';
            else
                cellsToUse = ~ROIinArtifact'  & ~offTargetRisk(holo,:);
            end
            popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if holo~=0
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
                
                distBins = [0:25:1000];
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
highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<0.1)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));

excludeInds = ismember(ensIndNumber,[]); %



ensemblesToUse = numSpikesEachEns > 75 &...
    numSpikesEachEns <110 &...
    lowRunInds &...
    ensStimScore > 0.5 &... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    ~excludeInds ;%&...  %
     %& numCellsEachEns>10 ;

indsSub = ensIndNumber(ensemblesToUse);
IndsUsed = unique(ensIndNumber(ensemblesToUse));

sum(ensemblesToUse)