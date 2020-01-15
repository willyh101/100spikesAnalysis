%% Load Experiments

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

numExps = numel(loadList);

clear All
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% Known Error Manual Catching
All(1).out.exp.stimParamsBackup = All(1).out.exp.stimParams;
All(1).out.exp.holoTargetsBackup = All(1).out.exp.holoTargets;
%%
All(1).out.exp.stimParams.Seq([3 4])=[];
All(1).out.exp.stimParams.numPulse([3 4])=[];
All(1).out.exp.stimParams.roi([3 4])=[];
All(1).out.exp.stimParams.Hz([2 3])=[];
All(1).out.exp.stimParams.numCells([2 3])=[];
All(1).out.exp.holoTargets([2 3])=[];



%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=6;
recWinSec = [1.25 2.5];

 for ind =1:numExps
     pTime =tic;
     fprintf(['Processing Experiment ' num2str(ind) '...']);
     
     All(ind).out.anal.numCells = size(All(ind).out.exp.zdfData,1);
     numCells(ind) = size(All(ind).out.exp.zdfData,1);
     
     if ~isfield(All(ind).out.info,'FR')
         All(ind).out.info.FR=FRDefault;
     end
     
     sz = size(All(ind).out.exp.zdfData);
     
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
     rdata = squeeze(mean(All(ind).out.exp.zdfData(:,winToUse,:),2));
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);
     bdata = squeeze(mean(All(ind).out.exp.zdfData(:,bwinToUse,:),2));
     
     All(ind).out.exp.rdData=rdata;
     All(ind).out.exp.bdata=bdata;

     sz2 = size(All(ind).out.vis.zdfData);
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz2(2)]) ;
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);

     rdata = squeeze(mean(All(ind).out.vis.zdfData(:,winToUse,:),2));
     bdata = squeeze(mean(All(ind).out.vis.zdfData(:,bwinToUse,:),2));
     
     All(ind).out.vis.rdata=rdata;
     All(ind).out.vis.bdata=bdata;
     
     temp = unique([All(ind).out.exp.holoTargets{:}]);
     temp(isnan(temp))=[];
     All(ind).out.anal.targets = temp;
     numUniqueTargets(ind) =numel(temp);
     
     fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

 end
%% Determine the OSI from the Vis section of each cell.

for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    uVisID = unique(All(ind).out.vis.visID);
    uVisID(uVisID==0)=[];
    
    oriCurve=[];
    for i=1:numel(uVisID)
        v= uVisID(i);
        trialsToUse = All(ind).out.vis.visID==v & All(ind).out.vis.lowMotionTrials;
        
        oriCurve(i,:)=mean(All(ind).out.vis.rdata(:,trialsToUse),2);
    end
    
    All(ind).out.anal.oriCurve = oriCurve;
    
    [maxOriVal maxOriIndex]= max(oriCurve);
    All(ind).out.anal.prefOri = maxOriIndex;
    
    % (Rpref ? Rorth)/(Rpref + Rorth)
    % visID 1 is always catch (I Hope... Will is more confident...ish)
    
    prefOri = maxOriIndex;
    orthoOri = prefOri-2;
    orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
    
    orthoOri2 = orthoOri+4;
    orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
    
    orthoOri = cat(1,orthoOri, orthoOri2);
    
    
    
    %     orthoOri(prefOri==1)=NaN;
    
    oriCurveBL = oriCurve - min(oriCurve);%oriCurve(1,:);
    
    OSI=[];
    for i=1:numel(prefOri)
        OSI(i) = (oriCurveBL(prefOri(i),i) - mean(oriCurveBL(orthoOri(:,i)',i)) ) / (oriCurveBL(prefOri(i),i)+ mean(oriCurveBL(orthoOri(:,i)',i)) );
        %     OSI = (oriCurveBL(prefOri) - oriCurveBL(orthoOri) ) ./ ( oriCurveBL(prefOri)+oriCurveBL(orthoOri) )
        OSI(prefOri==1)=nan;
    end
    
    All(ind).out.anal.OSI=OSI;
    
    
    pVisR=[];pVisT=[];
    for i=1:All(ind).out.anal.numCells
        trialsToUse = All(ind).out.vis.visID~=0 & All(ind).out.vis.lowMotionTrials;
        pVisR(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
        
        trialsToUse = All(ind).out.vis.visID~=0 & All(ind).out.vis.visID~=1 & All(ind).out.vis.lowMotionTrials;
        pVisT(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
    end
    
    All(ind).out.anal.pVisR = pVisR;
    All(ind).out.anal.pVisT = pVisT;
    
    alpha = 0.05;
    
    All(ind).out.anal.visPercent = sum(pVisR<alpha) / numel(pVisR);
    visPercent(ind) =  All(ind).out.anal.visPercent;

    meanOSI=[];ensembleOSI=[];ensembleOriCurve =[];ensemblePref=[];
    
    for i=1:numel(All(ind).out.exp.holoTargets)
        ht = All(ind).out.exp.holoTargets{i};
        ht(isnan(ht))=[];
        meanOSI(i)=nanmean(OSI(ht));
        
        ensOriCurve = mean(oriCurve(:,ht),2);
        ensembleOriCurve(i,:)= ensOriCurve;
        [maxOriVal maxOriIndex]= max(ensOriCurve);
        prefOri = maxOriIndex;
        
        orthoOri = prefOri-2;
        orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
        orthoOri2 = orthoOri+4;
        orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
        orthoOri = cat(1,orthoOri, orthoOri2);
        oriCurveBL = ensOriCurve - min(ensOriCurve);
        
        ensembleOSI(i) = (oriCurveBL(prefOri)- mean(oriCurveBL(orthoOri'))) / (oriCurveBL(prefOri)+ mean(oriCurveBL(orthoOri')));
        ensemblePref(i) = prefOri;
    end
    
    deg = [nan 0:45:315];
    All(ind).out.anal.ensembleOSI=ensembleOSI;
    All(ind).out.anal.meanOSI=meanOSI;
    All(ind).out.anal.ensembleOSI=ensembleOSI;
    All(ind).out.anal.ensembleOriCurve=ensembleOriCurve;
    All(ind).out.anal.ensemblePrefOri=ensemblePref;
    All(ind).out.anal.ensemblePrefDeg=deg(ensemblePref);
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% Pretty plots of OSI and tunings

clear allOSI ensOSI meanOSI ensNum roiNum
% OSI across all cells, all experiments
for i = 1:numel(All)
    allOSI{i} = All(i).out.anal.OSI(:);
    ensOSI{i} = All(i).out.anal.ensembleOSI(:);
    meanOSI{i} = All(i).out.anal.meanOSI(:);
    
    ensNum{i} = cellfun(@(x) sum(~isnan(x)),All(i).out.exp.holoTargets)'; %number of discovered Cells in ensemble
    roiNum{i} = cellfun(@(x) sum(~isnan(x)),All(i).out.exp.rois)'; %number of shot targets in ensemble
end

% unroll
allOSI = cell2mat(allOSI(:));
ensOSI = cell2mat(ensOSI(:));
meanOSI = cell2mat(meanOSI(:));
ensNum = cell2mat(ensNum(:));
roiNum = cell2mat(roiNum(:));

% plot
figure(1)
clf(1)
colors = {rgb('royalblue'), rgb('firebrick'), rgb('coral')};

hold on
h(1) = histogram(allOSI, 25);
h(2) = histogram(ensOSI, 25);
h(3) = histogram(meanOSI, 25);

for i=1:numel(h)
    h(i).Normalization = 'pdf';
    h(i).BinWidth = 0.05;
    h(i).FaceColor = colors{i};
    h(i).FaceAlpha = 0.44;
    kde(i) = fitdist(h(i).Data, 'kernel');
    p = plot(h(i).BinEdges, pdf(kde(i),h(i).BinEdges));
    p.LineWidth=2;
    p.Color= colors{i};
end

hold off

set(gcf(),'Name','OSI Across All Expts')
title('OSI Across All Expts')
xlabel('Orientation Selectivity Index')
ylabel('PDF')
legend('All Cells', 'Ensemble', 'Ensemble Mean')

%% now look at tuned ensembles
% get plots for the 2 different methods with higher bin count
f2 = figure(2);
clf(f2)
hold on
h(1) = histogram(ensOSI, 50);
%h(2) = histogram(meanOSI, 50);

% try using ensembleOSI > 0.3 for "tuned", this is arbitrary
OSIthreshold = 0.3;
for i = 1:numel(All)
    istuned = All(i).out.anal.ensembleOSI > OSIthreshold;
    tunedEnsembles = All(i).out.exp.holoTargets(istuned);
    untunedEnsembles = All(i).out.exp.holoTargets(~istuned);
    tunedEnsembleIdx = find(All(i).out.anal.ensembleOSI >= OSIthreshold);
    untunedEnsembleIdx = find(All(i).out.anal.ensembleOSI < OSIthreshold);
    All(i).out.anal.tunedEnsembles = tunedEnsembles;
    All(i).out.anal.untunedEnsembles = untunedEnsembles;
    All(i).out.anal.tunedEnsembleIdx = tunedEnsembleIdx;
    All(i).out.anal.untunedEnsembleIdx = untunedEnsembleIdx;
end

% plot the OSIs of tuned vs untuned ensembles
for i = 1:numel(All)
    allOSIT{i} = All(i).out.anal.OSI(:);
    ensOSIT{i} = All(i).out.anal.ensembleOSI(:);
    meanOSIT{i} = All(i).out.anal.meanOSI(:);
end

% unroll
allOSIT = cell2mat(allOSIT(:));
ensOSIT = cell2mat(ensOSIT(:));
meanOSIT = cell2mat(meanOSIT(:));

f3 = subplots(1,2,1)

% and thier preferred oris

    
    
%% Get the number of spikes in each stimulus

clear numSpikesEachStim numCellsEachEns
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
end
numSpikesEachStim=cell2mat(numSpikesEachStim(:)');
numSpikesEachEns = numSpikesEachStim;
numSpikesEachEns(numSpikesEachStim==0)=[];

numCellsEachEns=cell2mat(numCellsEachEns(:)');
    
%% Make all dataPlots into matrixes of mean responses


clear popResponse
for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);

    trialsToUse = All(ind).out.exp.lowMotionTrials;

    clear respMat baseMat %Order stims,vis,cells
    for i=1:numel(unique(All(ind).out.exp.stimID))
        us = unique(All(ind).out.exp.stimID);
        s = us(i);

        for k= 1 : numel(unique(All(ind).out.exp.visID))
            vs = unique(All(ind).out.exp.visID);
            v = vs(k);

            respMat(i,k,:) = mean(All(ind).out.exp.rdData(:,...
                trialsToUse & All(ind).out.exp.stimID ==s &...
                All(ind).out.exp.visID ==v), 2) ;
            baseMat(i,k,:) = mean(All(ind).out.exp.bdata(:,...
                trialsToUse & All(ind).out.exp.stimID ==s &...
                All(ind).out.exp.visID ==v), 2) ;
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

    thisPlaneTolerance = 10;10; %in pixels
    onePlaneTolerance = 20;20;

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

    %%Get Pop Responses
    %         v=1; %best bet for no vis stim.
    clear popResp
    for v = 1:numel(vs)
        for i= 1:numel(All(ind).out.exp.stimParams.Seq)
            %             try
            holo =All(ind).out.exp.stimParams.Seq(i) ;% roi{i}{1};
            %             catch
            %                 holo =All(ind).out.exp.stimParams.roi{i};
            %             end

            if i==1;
                cellsToUse = ~ROIinArtifact';
            else
                cellsToUse = ~ROIinArtifact' & ~offTargetRisk(holo,:);
            end
            popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
        end
    end

    popResponse{ind} = popResp(:,1);
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

popResponse = cell2mat(popResponse(:));
popResponseEns=popResponse;
popResponseEns(numSpikesEachStim==0)=[];
    
%% Plot
figure(37);
ensemblesToUse = numSpikesEachEns>75 ;%& numCellsEachEns<10 ;
scatter(meanOSI(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')

xlabel('Ensemble OSI')
ylabel('Population Mean Response')
    
%% what do the mean/normalized tuning look like for low, middle, and high OSI values?
clear oriShifted lowOSIidx midOSIidx highOSIidx lowOSIcurve midOSIcurve highOSIcurve alignedOris
% set bounds for OSI
low = 0.4;
high = 0.6;

for i = 1:numel(All)
    
    lowOSIidx{i} = find(All(i).out.anal.ensembleOSI <= low);
    midOSIidx{i} = find(All(i).out.anal.ensembleOSI > low & All(i).out.anal.ensembleOSI < high);
    highOSIidx{i} = find(All(i).out.anal.ensembleOSI >= high);
    
    ensembleOriCurve = All(i).out.anal.ensembleOriCurve;
    ensemblePref = All(i).out.anal.ensemblePrefOri;
    
    % peak align to the 3rd position
    % should this be normalized somehow?
    for j = 1:size(ensembleOriCurve,1)
        oriShifted(j,:) = circshift(ensembleOriCurve(j,:),-ensemblePref(j)+3);
    end
    
    alignedOris{i} = oriShifted;
    
    lowOSIcurve{i} = oriShifted(lowOSIidx{i},:);
    midOSIcurve{i} = oriShifted(midOSIidx{i},:);
    highOSIcurve{i} = oriShifted(highOSIidx{i},:);
    
    All(i).out.anal.lowOSIcurve = oriShifted(lowOSIidx{i},:);
    All(i).out.anal.midOSIcurve = oriShifted(midOSIidx{i},:);
    All(i).out.anal.highOSIcurve = oriShifted(highOSIidx{i},:);

end

% unroll and computer errors
lowOSIcurveAll = cell2mat(lowOSIcurve(:));
err1 = std(lowOSIcurveAll)/sqrt(size(lowOSIcurveAll,1));
midOSIcurveAll = cell2mat(midOSIcurve(:));
err2 = std(midOSIcurveAll)/sqrt(size(midOSIcurveAll,1));
highOSIcurveAll = cell2mat(highOSIcurve(:));
err3 = std(highOSIcurveAll)/sqrt(size(highOSIcurveAll,1));

% plot them
f4 = figure(4);
clf(4)
hold on
errorbar(nanmean(lowOSIcurveAll,1), err1, 'linewidth',2);
errorbar(nanmean(midOSIcurveAll,1), err2, 'linewidth', 2);
errorbar(nanmean(highOSIcurveAll,1), err3, 'linewidth', 2);
hold off

title('Mean OSI Curves')
ylabel('Mean Response')
xlabel('Ori (preferred centered at 3)')
legend('Low OSIs', 'Mid OSIs', 'High OSIs')
set(gcf(),'Name','Mean OSI Curves')






