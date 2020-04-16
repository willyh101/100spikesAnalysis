function [All, outVars] = meanMatrixVisandCorr(All,opts,outVars);
%% Make all dataPlots into matrixes of mean responses
 %%Determine Vis Responsive and Process Correlation
visAlpha =  opts.visAlpha;% 0.05;
 
 %oftarget risk params
 thisPlaneTolerance = opts.thisPlaneTolerance;%11.25;%7.5;%1FWHM%10; %in um;% pixels
 onePlaneTolerance = opts.onePlaneTolerance;% 22.5;%15;%2FWHM %20;
 
numExps = numel(All);

numSpikesEachStim = outVars.numSpikesEachStim;

clear popResponse pVisR pVisT
ensIndNumber=[];
ensHNumber=[];
for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);

    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;
    
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);

    clear respMat baseMat stdMat numMat%Order stims,vis,cells
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
                stdMat(i,v,:) = std(All(ind).out.exp.rdData(:,...
                    trialsToUse & All(ind).out.exp.stimID ==s &...
                    All(ind).out.exp.visID ==v),[], 2) ;
                numMat(i,v,:) = sum(trialsToUse & All(ind).out.exp.stimID ==s &...
                    All(ind).out.exp.visID ==v) ;
            end
        end
    end

    All(ind).out.anal.respMat = respMat;
    All(ind).out.anal.baseMat = baseMat;
    All(ind).out.anal.stdMat = stdMat;
    All(ind).out.anal.numMat = numMat;


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
    All(ind).out.anal.pVisR = pVisR;
    
     % hack for expt 11
    % use the vis stim expt
%     if ind==11
%         pVisR=[];
%         for i=1:All(ind).out.anal.numCells     
%             trialsToUse = All(ind).out.vis.visID~=0 &...
%                 All(ind).out.vis.lowMotionTrials &...
%                 All(ind).out.vis.lowRunTrials;
%             pVisR(i) = anova1(All(ind).out.vis.rdData(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
%         end
%     end
%     All(ind).out.anal.pVisR = pVisR;
     
    
    All(ind).out.anal.visPercent = sum(pVisR<visAlpha) / numel(pVisR);
    visPercent(ind) =  All(ind).out.anal.visPercent;

    %%Get Pop Responses
    %         v=1; %best bet for no vis stim.
    vs(vs==0)=[];
    clear popResp popRespDist popRespDistNumCells popRespDistSubtracted  popRespDistVisNumCells popRespDistSubVis popRespDistVis
    clear minDistbyHolo geoDistbyHolo meanDistbyHolo harmDistbyHolo
    
    
    distBins = opts.distBins; %[0:25:1000];
    numDist =numel(distBins)-1;
    numStims = numel(us); %numel(All(ind).out.exp.stimParams.Seq); %Caution I made this change but i'm not sure about it -Ian
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
            popResp(i,v) = mean(squeeze(respMat(i,v ,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if i~=1
                Tg=All(ind).out.exp.rois{holo};
                dists = StimDistance(Tg,:);
                
                minDist = min(dists,[],1);
                try
                    geoDist = geo_mean(dists,1);
                catch
                    geoDist = geomean(dists,1);
                end
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
    ensHNumber = [ensHNumber 1:numel(popResp(:,1))];

    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

popResponse = cell2mat(popResponse(:));
popResponseEns=popResponse;
popResponseEns(numSpikesEachStim==0)=[];

ensIndNumber(numSpikesEachStim==0)=[];
ensHNumber(numSpikesEachStim==0)=[];

noStimPopResp = popResponse(numSpikesEachStim==0);

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<0.05)); %remove low vis responsive experiments

outVars.ensIndNumber            = ensIndNumber; %the ind that each holo came from
outVars.ensHNumber            = ensHNumber; %the holo ID number index of uniquestims

outVars.visPercent              = visPercent;
outVars.popResponse             = popResponse;
outVars.popResponseDist         = popResponseDist;
outVars.popResponseNumCells     = popResponseNumCells;
outVars.popResponseDistVis      = popResponseDistVis;
outVars.popResponseEns          = popResponseEns;
outVars.ensIndNumber            = ensIndNumber;
outVars.noStimPopResp           = noStimPopResp;
outVars.highVisPerentInd        = highVisPercentInd;

% outVars.Others = {popResponseAll;...
%     popResponseAllDist;...
%     popResponseAllDistSub;...
%     popResponseAllNumCells;...
%     popResponseAllDistSubVis;...
%     popResponseAllDistVis;...
%     popResponseAllDistSubVisNC};

outVars.popResponseAll          = popResponseAll;
outVars.popResponseAllDist      = popResponseAllDist;
outVars.popResponseAllDistSub   = popResponseAllDistSub;
outVars.popResponseAllNumCells  = popResponseAllNumCells;
outVars.popResponseAllDistSubVis = popResponseAllDistSubVis;
outVars.popResponseAllDistVis   = popResponseAllDistVis;
outVars.popResponseAllDistSubVisNC = popResponseAllDistSubVisNC;



