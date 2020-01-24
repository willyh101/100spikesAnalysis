%% Load Experiments

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

% loadPath = 'Z:\ioldenburg\outputdata1'
%%
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
     
     %ensure has a visID
     if ~isfield(All(ind).out.exp,'visID')
         All(ind).out.exp.visID = ones(size(All(ind).out.exp.stimID));
         disp(['Added visID to Exp ' num2str(ind)]);
     end
     
     if numel(All(ind).out.vis.visID) ~= numel(All(ind).out.vis.lowMotionTrials)
         All(ind).out.vis.lowMotionTrials(end+1:numel(All(ind).out.vis.visID))= 0 ;
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
     
     %infer a visStart 
%      All(7).out.exp.outputsInfo.OutputPatterns{2}
     
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
    
    visAlpha = 0.05;
    
    All(ind).out.anal.visPercent = sum(pVisR<visAlpha) / numel(pVisR);
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

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<0.05)); %remove low vis responsive experiments


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
h(1) = histogram(meanOSI, 50);
h(2) = histogram(ensOSI, 50);


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
clear allOSIT ensOSIT meanOSIT
for i = 1:numel(All)
    allOSIT{i} = All(i).out.anal.OSI(:);
    ensOSIT{i} = All(i).out.anal.ensembleOSI(:);
    meanOSIT{i} = All(i).out.anal.meanOSI(:);
end

% unroll
allOSIT = cell2mat(allOSIT(:));
ensOSIT = cell2mat(ensOSIT(:));
meanOSIT = cell2mat(meanOSIT(:));

% f3 = subplots(1,2,1)

% and thier preferred oris

legend('Mean OSI', 'Ensemble OSI')
    
%% Get the number of spikes in each stimulus

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
%% Make all dataPlots into matrixes of mean responses


clear popResponse pVisR pVisT
ensIndNumber=[];
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
    holoTargets = All(ind).out.exp.holoTargets;

    thisPlaneTolerance = 10;10; %in pixels
    onePlaneTolerance = 15;20;

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
    pVisR = All(ind).out.anal.pVisR;
    pVisT = All(ind).out.anal.pVisT;

    %%Get Pop Responses
    %         v=1; %best bet for no vis stim.
    clear popResp popRespDist popRespDistNumCells 
    clear minDistbyHolo geoDistbyHolo meanDistbyHolo harmDistbyHolo
    for v = 1:numel(vs)
        for i= 1:numel(All(ind).out.exp.stimParams.Seq)
            %             try
%             holo =All(ind).out.exp.stimParams.Seq(i) ;% roi{i}{1};
            holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
            %             catch
            %                 holo =All(ind).out.exp.stimParams.roi{i};
            %             end

            if i==1;
                cellsToUse = ~ROIinArtifact' & pVisR<0.05;
            else
                cellsToUse = ~ROIinArtifact' & pVisR<0.05 & ~offTargetRisk(holo,:);
            end
            popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if i~=1
                Tg=All(ind).out.exp.rois{holo};
                dists = StimDistance(Tg,:);
                minDist = min(dists,[],1);
                geoDist = geo_mean(dists,1); 
                meanDist = mean(dists,1);
                harmDist = harmmean(dists,1);
                
                minDistbyHolo(i,:) = minDist;
                geoDistbyHolo(i,:) = geoDist;
                meanDistbyHolo(i,:) = meanDist;
                harmDistbyHolo(i,:) = harmDist;
                
                distToUse = minDist; % CHANGE THIS (when you want to change whats being analyzed)
                
                distBins = [0:25:500];
                for d = 1:numel(distBins)-1
                    cellsToUse = ~ROIinArtifact' &...
                        pVisR<10.05 &... %% IS VIS RESPONSIVE?
                        ~offTargetRisk(holo,:) &...
                        distToUse > distBins(d) &...
                        distToUse <= distBins(d+1) ;
                    popRespDist(i,v,d) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                    popRespDistNumCells(i,v,d) = sum(cellsToUse);
                end
            end
        
        
        end
    end
    
    VisCondToUse = 1; %1 is no vis
    if VisCondToUse > size(popResp,2) 
        popResponse{ind} = single(nan(size(popResp(:,1))));
        popResponseDist{ind} = single(nan(size(squeeze(popRespDist(:,1,:)))));
        popResponseNumCells{ind} = single(nan(size(squeeze(popRespDistNumCells(:,1,:)))));
    else
        popResponse{ind} = popResp(:,VisCondToUse);
        popResponseDist{ind} = squeeze(popRespDist(:,VisCondToUse,:));
        popResponseNumCells{ind} = squeeze(popRespDistNumCells(:,VisCondToUse,:));
    end
    popResponseAll{ind} = popResp;
    popResponseAllNumCells{ind} = popRespDistNumCells;
    
    ensIndNumber = [ensIndNumber ones(size(popResp(:,1)'))*ind];
        All(ind).out.anal.minDistbyHolo = minDistbyHolo;
        All(ind).out.anal.geoDistbyHolo = geoDistbyHolo;
        All(ind).out.anal.meanDistbyHolo = meanDistbyHolo;

    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

popResponse = cell2mat(popResponse(:));
popResponseEns=popResponse;
popResponseEns(numSpikesEachStim==0)=[];

ensIndNumber(numSpikesEachStim==0)=[];

noStimPopResp = popResponse(numSpikesEachStim==0);

    
%% Plot
f3 = figure(3);
clf(3)


ensemblesToUse = numSpikesEachEns > 75 & numSpikesEachEns <125 & highVisPercentInd;% & ensIndNumber==15; %& numCellsEachEns>10 ;
%scatter(meanOSI(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
scatter(ensOSI(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')

% p = polyfit(ensOSI(ensemblesToUse),popResponseEns(ensemblesToUse),1);
% f = polyval(p, ensOSI(ensemblesToUse));
% hold on
% plot(ensOSI(ensemblesToUse), f)
% hold off
xlabel('Ensemble OSI')
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

%% fit ensembles of different sizes
clear f p ens2plt fits
f5 = figure(5);
clf(f5)
numEns = numel(unique(numCellsEachEns(ensemblesToUse)));

uniqueEns = unique(numCellsEachEns(ensemblesToUse));

for i=1:numEns
    ens2plot = find(numCellsEachEns(ensemblesToUse)==uniqueEns(i));
    p = polyfit(ensOSI(ens2plot),popResponseEns(ens2plot),1);
    f = polyval(p, ensOSI(ens2plot));
    %fits(i,:) = f;
    subplot(1,numEns,i)
    plt = plot(ensOSI(ens2plot),popResponseEns(ens2plot), '.', 'MarkerSize',12);
    hold on
    fline = plot(ensOSI(ens2plot), f, 'LineWidth', 1);
    xlabel('OSI')
    ylabel('Pop Response')
    title(['Ensembles of size ' num2str(uniqueEns(i))])
end

linkaxes

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

ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))


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

%% what do the mean/normalized tuning look like for low, middle, and high OSI values?
clear oriShifted lowOSIidx midOSIidx highOSIidx lowOSIcurve midOSIcurve highOSIcurve alignedOris
% set bounds for OSI
low = 0.2;
high = 0.8;

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



%% Plot Pop Response by Distance
popDistAll = cell2mat(popResponseDist');
popDist = popDistAll(numSpikesEachStim~=0,:);

popNumCells = cell2mat(popResponseNumCells');
popNCells = popNumCells(numSpikesEachStim~=0,:);


ensSizes = unique(numCellsEachEns(ensemblesToUse))   ;


% colorList = {rgb('DarkBlue') rgb('steelblue') rgb('gold')};
colorList = colormap(viridis(numel(ensSizes)));
colorList = num2cell(colorList,2);

figure(8);clf
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDist(ensemblesToUse & numCellsEachEns==ensSizes(i) & highVisPercentInd,:);
meanDat = nanmean(dat)
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
% xlim([0 400])
legend('Small', 'Medium', 'Big')
legend(string(ensSizes))

%% Ensemble stims and vis things

% for i=1:numExps
%     
% 
% 



%% Correlation Analyisis Determine Correlation Coefficients

for ind = 1:numExps
      pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
%     us = unique(All(ind).out.vis.stimID);

    %Spont Corr - correlation coefficient on time series from no stim
    %period
    trialsToUse = All(ind).out.vis.lowMotionTrials &...
        All(ind).out.vis.visID == 1;
    unrollData = All(ind).out.vis.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    [SpontCorr SpCoP] = corr(unrollData');
    
    %AllCorr - the correlation coef on all time series
    trialsToUse = All(ind).out.vis.lowMotionTrials ;
    unrollData = All(ind).out.vis.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    [AllCorr AlCoP] = corr(unrollData');
    
    %All corr mean - correlation coef of response (not time series)
    trialsToUse = All(ind).out.vis.lowMotionTrials;
    unrollData = All(ind).out.vis.rdata(:,trialsToUse);
    sz = size(unrollData);
    
    [AllMCorr AmCoP] = corr(unrollData');
    
    %noise corr - correlation coef of residual trial response (not time
    %series) i.e. trial response - mean response for that condition
     trialsToUse = All(ind).out.vis.lowMotionTrials ;
    vs = unique(All(ind).out.vis.visID);
    vs(vs==0)=[];
    unrollData = [];
    meanResps = [];
    for k = 1:numel(vs)
        v = vs(k);
        trialsToUseThis = trialsToUse & All(ind).out.vis.visID==v;
        
        dataPart = All(ind).out.vis.rdata(:,trialsToUseThis);
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

%% plot Pop Response by Correlation
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
    if all(All(ind).out.exp.visID==0)
        All(ind).out.exp.visID = ones(size(All(ind).out.exp.visID));
    else
        vs(vs==0)=[];
    end
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    clear popRespCorr popRespCorrSub minDistbyHolo cellsToUse
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
                meanCorr = mean(distCorr,1);
%                 geoCorr = geomean(distCorr,1);
                minCorr = min(distCorr,[],1);
                maxCorr = max(distCorr,[],1);
                
                corrHolo = meanCorr;
              
                
                distBins = linspace(-1,1,80);
                for d = 1:numel(distBins)-1
                    cellsToUse = ~ROIinArtifact' &...
                        ~offTargetRisk(holo,:) &...
                        corrHolo > distBins(d) &...
                        corrHolo <= distBins(d+1) ;
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
 
%  All(ind).out.anal.minDistbyHolo = minDistbyHolo;

end


% 
% popResponseCorr = cell2mat(popResponseCorr(:));
% popResponseCorr(numSpikesEachStim==0)=[];

popRespCorr;


% temp = cellfun(@(x) squeeze(x(:,1,:)),popResponseCorr,'uniformoutput',0) ;

temp = cellfun(@(x) reshape(x(:,1,:),size(x(:,1,:),1),size(x(:,1,:),3)),popResponseCorr,'uniformoutput',0) ;
% temp = cellfun(@(x) reshape(x(:,end,:),size(x(:,end,:),1),size(x(:,end,:),3)),popResponseCorrSub,'uniformoutput',0) ;

% temp = cellfun(@(x) squeeze(x(:,end,:)),popResponseCorr,'uniformoutput',0) ;

% temp = cellfun(@(x) squeeze(x(:,end,:)),popResponseCorrSub,'uniformoutput',0) ;
% temp = cellfun(@(x) squeeze(x(:,1,:)),popResponseCorrSub,'uniformoutput',0) ;

% temp = cellfun(@(x) squeeze(x(:,round(end/2),:)),popResponseCorr,'uniformoutput',0) ;

EnsCorR = cat(1,temp{:});

figure(15);clf
% subplot(1,3,1)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse & highVisPercentInd);
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
% xlim([-0.4 0.4])
legend('small', 'medium', 'large')
r = refline(0);
r.LineStyle = ':';
r.Color = rgb('grey');
r.LineWidth = 2;

ylabel('Pop Resp to HoloStim')
xlabel('Responder to Ensemble Correlation')



%%

%% Calculate L1 and L2 
contrastsToView = [6 3 2 1.5 1.25 1] ;%I know its weird i just wanted to be able to catch times that we were using different contrasts, will work out to 1:6 if there are 6 contrasts; 1:6;

clear EnsL1 EnsL2
c=0;
for ind =1:numExps
    
  for h= 1:numel(All(ind).out.exp.stimParams.Seq)-1
            holo = All(ind).out.exp.stimParams.roi{h+1}; % only cycle through holos
    
            divider =1;
            maxV = max(All(ind).out.exp.visID);
            v = max(round(maxV/divider),1);
            
            trialsToUse = All(ind).out.exp.lowMotionTrials & All(ind).out.exp.visID==v;
%             cellsToUse =  ~All(ind).out.anal.ROIinArtifact' & All(ind).out.anal.offTargetRisk(holo,:);
                        cellsToUse =  ~All(ind).out.anal.ROIinArtifact' & ~any(All(ind).out.anal.offTargetRisk(:,:));

                        
            us = unique(All(ind).out.exp.stimID); 
            
            testData = All(ind).out.exp.rdData(cellsToUse,trialsToUse & All(ind).out.exp.stimID == us(h+1));
            ExpectedData = All(ind).out.exp.rdData(cellsToUse,trialsToUse & All(ind).out.exp.stimID == us(1));
            [L1 L2 L3] =  calcL1L2(testData,ExpectedData);
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
clear data names
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{i} = EnsL2(ens2plot)./EnsL1(ens2plot);
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
ylabel('L2/L1')
ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

% pValEnselbeSize = anovan(popResponseEns(ensemblesToUse),numCellsEachEns(ensemblesToUse)','display','off')
% 
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))
%% L1 L2 by dist
                distBins = [0:25:1000];
                
                numEns = numel(uniqueEns);
                
                clear EnsL1 EnsL2
c=0;
for ind =1:numExps
    
  for h= 1:numel(All(ind).out.exp.stimParams.Seq)-1
            holo = All(ind).out.exp.stimParams.roi{h+1}; % only cycle through holos
%             divider = 1; %6 is no vis, 1 is max vis
%             x = max(All(ind).out.exp.visID);
%             v = max(round(size(x,2)/divider),1);
            
            trialsToUse = All(ind).out.exp.lowMotionTrials & All(ind).out.exp.visID==1;
            cellsToUse =  ~All(ind).out.anal.ROIinArtifact' & ~All(ind).out.anal.offTargetRisk(holo,:);
            
            us = unique(All(ind).out.exp.stimID); 
            
            c=c+1;

            distByHoloToUse =  All(ind).out.anal.geoDistbyHolo;
            
            for d = 1:numel(distBins)-1
                D = distBins(d+1);
                cellsToUseDist = cellsToUse &...
                   distByHoloToUse(h+1,:) <=distBins(d+1) &...
                    distByHoloToUse(h+1,:) >distBins(d) ;
            
            testData = All(ind).out.exp.rdData(cellsToUseDist,trialsToUse & All(ind).out.exp.stimID == us(h+1));
            ExpectedData = All(ind).out.exp.rdData(cellsToUseDist,trialsToUse & All(ind).out.exp.stimID == us(1));
            [L1 L2] =  calcL1L2(testData,ExpectedData);
            L1= L1/size(testData,1);
            L2 = L2/sqrt(size(testData,1));
            
            EnsL1(c,d) = L1;
            EnsL2(c,d) = L2;   
            
            end
  end
    
end



figure(111);clf
subplot(1,3,1)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse & highVisPercentInd);
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
% xlim([0 250])
legend('small', 'medium', 'large')
title('L1')
xlabel('Distance')

subplot(1,3,2)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse & highVisPercentInd);
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
% xlim([0 250])
legend('small', 'medium', 'large')
title('L2')
xlabel('Distance')

subplot(1,3,3)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse & highVisPercentInd);
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
% xlim([0 250])
legend('small', 'medium', 'large')
title('L2/L1 Ratio')
xlabel('Distance')
