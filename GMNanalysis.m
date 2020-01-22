%% Load Experiments

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

% 'U:\ioldenburg\outputdata1'
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

%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=6;
recWinRange = [1 2];% %from vis Start [1.25 2.5];

 for ind =1:numExps
     pTime =tic;
     fprintf(['Processing Experiment ' num2str(ind) '...']);
     
     All(ind).out.anal.numCells = size(All(ind).out.exp.zdfData,1);
     numCells(ind) = size(All(ind).out.exp.zdfData,1);
     
     if ~isfield(All(ind).out.info,'FR')
         All(ind).out.info.FR=FRDefault;
     end
     
     sz = size(All(ind).out.exp.zdfData);
         
     
     if ~isfield(All(ind).out.exp,'visStart')
         All(ind).out.exp.visStart = 0.5;
         disp(['Added visStart to Exp' num2str(ind)]);
     end
     
     recWinSec = recWinRange + All(ind).out.exp.visStart; %recording window relative to when vis start
     
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
     rdata = squeeze(mean(All(ind).out.exp.zdfData(:,winToUse,:),2));
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);
     bdata = squeeze(mean(All(ind).out.exp.zdfData(:,bwinToUse,:),2));
     
     All(ind).out.exp.rdData=rdata;
     All(ind).out.exp.bdata=bdata;

    
     temp = unique([All(ind).out.exp.holoTargets{:}]);
     temp(isnan(temp))=[];
     All(ind).out.anal.targets = temp;
     numUniqueTargets(ind) =numel(temp);
     
     %ensure has a visID
     if ~isfield(All(ind).out.exp,'visID')
         All(ind).out.exp.visID = ones(size(All(ind).out.exp.stimID));
         disp(['Added visID to Exp ' num2str(ind)]);
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
     
     
     fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

 end
 

 
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
 %%Determine Vis Responsive and Process Correlation
 

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
%     pVisR = All(ind).out.anal.pVisR;
%     pVisT = All(ind).out.anal.pVisT;


    %ID tuned Cells, should comparing no contrast to with contrast
    pVisR=[];%pVisT=[];
    for i=1:All(ind).out.anal.numCells
        trialsToUse = All(ind).out.exp.visID~=0 & All(ind).out.exp.lowMotionTrials & All(ind).out.exp.stimID==min(All(ind).out.exp.stimID);
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
    clear minDistbyHolo
    for v = 1:numel(vs)
        for i= 1:numel(All(ind).out.exp.stimParams.Seq)
            %             try
%             holo =All(ind).out.exp.stimParams.Seq(i) ;% roi{i}{1};
            holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
            %             catch
            %                 holo =All(ind).out.exp.stimParams.roi{i};
            %             end

            if i==1;
                cellsToUse = ~ROIinArtifact';
            else
                cellsToUse = ~ROIinArtifact'  & ~offTargetRisk(holo,:);
            end
            popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if i~=1
                Tg=All(ind).out.exp.rois{holo};
                dists = StimDistance(Tg,:);
                minDist = min(dists);
                
                minDistbyHolo(i,:) = minDist;
                distBins = [0:25:500];
                for d = 1:numel(distBins)-1
                    cellsToUse = ~ROIinArtifact' &...
                        ~offTargetRisk(holo,:) &...
                        minDist > distBins(d) &...
                        minDist <= distBins(d+1) ;
                    popRespDist(i,v,d) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                    popRespDistNumCells(i,v,d) = sum(cellsToUse);
                    popRespDistSubtracted(i,v,d) = popRespDist(i,v,d) - (mean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse))));
                    cellsToUse = cellsToUse & pVisR<visAlpha;
                    tempResp = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                    popRespDistSubVis(i,v,d) = tempResp - (mean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse))));
                    popRespDistVis(i,v,d) = tempResp;
                    popRespDistVisNumCells(i,v,d) = sum(cellsToUse);

                  
                end

            end
        
        
        end
    end
    
    All(ind).out.anal.minDistbyHolo = minDistbyHolo;
    
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
    popResponseAll{ind} = popResp;
    popResponseAllDist{ind} = popRespDist;
    popResponseAllDistSub{ind} = popRespDistSubtracted;
    popResponseAllNumCells{ind} = popRespDistNumCells;
    popResponseAllDistSubVis{ind} = popRespDistSubVis;
    popResponseAllDistVis{ind} = popRespDistVis;
    popResponseAllDistSubVisNC{ind} = popRespDistVisNumCells;
    
    ensIndNumber = [ensIndNumber ones(size(popResp(:,1)'))*ind];
    

    
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


ensemblesToUse = numSpikesEachEns > 75 & numSpikesEachEns <110;% & ensIndNumber==15; %& numCellsEachEns>10 ;
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


%% Plot Pop Response by Distance
popDistAll = cell2mat(popResponseDist');
popDist = popDistAll(numSpikesEachStim~=0,:);

popNumCells = cell2mat(popResponseNumCells');
popNCells = popNumCells(numSpikesEachStim~=0,:);


ensSizes = unique(numCellsEachEns(ensemblesToUse))   ;


colorList = {rgb('DarkBlue') rgb('steelblue') rgb('gold')};

figure(9);clf
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDist(ensemblesToUse & numCellsEachEns==ensSizes(i),:);
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
xlim([0 250])
legend('Small', 'Medium', 'Big')

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

subplot(1,2,1)


ensSizes = unique(numCellsEachEns(ensemblesToUse))   ;


colorList = {rgb('DarkBlue') rgb('steelblue') rgb('gold')};

for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatNoVis(ensemblesToUse & numCellsEachEns==ensSizes(i),:);
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
xlim([0 250])
legend('Small', 'Medium', 'Big')
title('Holographic induced changes 0 Contrast')

subplot(1,2,2)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatMaxVisSubtracted(ensemblesToUse & numCellsEachEns==ensSizes(i),:);
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
xlim([0 250])
legend('Small', 'Medium', 'Big')
title('Holographic induced changes Max Contrast')


figure(12);clf
contrastsToView = [6 3 2 1.5 1.25 1] ;%I know its weird i just wanted to be able to catch times that we were using different contrasts, will work out to 1:6 if there are 6 contrasts; 1:6;
for c=1:numel(contrastsToView)
ax(c) = subplot(1,numel(contrastsToView),c);
divider = contrastsToView(c);
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistSub,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDist,'uniformoutput',0)');
popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistSubVis,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistVis,'uniformoutput',0)');

for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatToPlot(ensemblesToUse & numCellsEachEns==ensSizes(i),:);
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
xlim([0 250])
legend('Small', 'Medium', 'Big')
title(['Contrast ' num2str(c) ] );



end

linkaxes([ax(:)])


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
    data{i} = EnsL1(ens2plot)./EnsL2(ens2plot);
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

% pValEnselbeSize = anovan(popResponseEns(ensemblesToUse),numCellsEachEns(ensemblesToUse)','display','off')
% 
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
% ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))
%% L1 L2 by dist
All(ind).out.anal.minDistbyHolo;
                distBins = [0:25:500];
                
                clear EnsL1 EnsL2
c=0;
for ind =1:numExps
    
  for h= 1:numel(All(ind).out.exp.stimParams.Seq)-1
            holo = All(ind).out.exp.stimParams.roi{h+1}; % only cycle through holos
    divider = 6;
            x = max(All(ind).out.exp.visID);
            v = max(round(size(x,2)/divider),1);
            
            trialsToUse = All(ind).out.exp.lowMotionTrials & All(ind).out.exp.visID==1;
            cellsToUse =  ~All(ind).out.anal.ROIinArtifact' & ~All(ind).out.anal.offTargetRisk(holo,:);
            
            us = unique(All(ind).out.exp.stimID); 
            
            c=c+1;

            for d = 1:numel(distBins)-1
                D = distBins(d+1);
                cellsToUseDist = cellsToUse &...
                    All(ind).out.anal.minDistbyHolo(h+1,:) <=distBins(d+1) &...
                    All(ind).out.anal.minDistbyHolo(h+1,:) >distBins(d) ;
            
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

figure(11);clf
subplot(1,3,1)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
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
xlim([0 250])
legend('small', 'medium', 'large')
subplot(1,3,2)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
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
xlim([0 250])
legend('small', 'medium', 'large')
subplot(1,3,3)
hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data = EnsL1(ens2plot,:)./EnsL2(ens2plot,:);
    
    
    e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data)./sqrt(sum(~isnan(data))));
%         e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data));

    e.Color = colorList{i};
    e.LineWidth = 2;
%     names{i} = string(uniqueEns(i));
%     avg(i) = mean(popResponseEns(ens2plot));
%     err(i) = sem(popResponseEns(ens2plot));
%     ns(i) = numel(popResponseEns(ens2plot));
end
xlim([0 250])
legend('small', 'medium', 'large')


