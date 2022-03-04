%% stand alone plot of some example Data

meanEnsResp = mean(outVars.popResponseEns(ensemblesToUse));
stdEnsResp = std(outVars.popResponseEns(ensemblesToUse));

eligibleEns = outVars.popResponseEns<meanEnsResp & outVars.popResponseEns<meanEnsResp>(meanEnsResp-stdEnsResp) & ensemblesToUse';
eligibleEns =  ensemblesToUse';

egL = find(eligibleEns);
for i =1:sum(eligibleEns)
    %Nice ones 220302: 594 768 771 822 828 829 836 837 838 845 866 879 1005
    %1027 1030 1037 1053 1056 1058
    
    %Nice 220303: 398 399 594 768 822 829 837 838 839 845 866 877 879 885
    %929 983
    %229 230 232 537 560 594 7 68 769
    ens = egL(i); %953; 
    thisEnsScore = outVars.popResponseEns(ens) ;
    
    
    ind = outVars.ensIndNumber(ens);
    h = outVars.ensHNumber(ens);
    % sID = All(ind).out.exp.stimParams.Seq(h)
    sID = h; %does this hold true if the first isn't a no hologram?
    
    
    rawDat = All(ind).out.exp.dfData;
    rawDatForCellVal = All(ind).out.exp.zdfData;
    
    us = unique(All(ind).out.exp.stimID);
    
    trialsToUse1 = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(sID) & ...
        All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);
    
    dataToUse = rawDat(:,:,trialsToUse1); %order cells, frames, trials
    win = All(ind).out.anal.recWinUsed;
    bwin = All(ind).out.anal.bwinToUse;
    
    cellVal = mean(mean(rawDatForCellVal(:,win(1):win(2),trialsToUse1),2),3)';
    bVal = mean(mean(rawDatForCellVal(:,bwin(1):bwin(2),trialsToUse1),2),3)';
    cellVal = cellVal-bVal;
    
    trialsToUse2 = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) & ...
        All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);
    
    
    noStimDataToUse = rawDat(:,:,trialsToUse2);
    
    
    trialAverageData = mean(dataToUse,3);
    
    offH = All(ind).out.exp.stimParams.roi{h};
    thisOffTargetRisk = logical(All(ind).out.anal.offTargetRisk(offH,:));
    
    eligibleCells  = ~thisOffTargetRisk...
        & ~All(ind).out.anal.ROIinArtifact';
    %     & cellVal<mean(cellVal) ...
    %     & cellVal>mean(cellVal)-std(cellVal) ...
    %     & cellVal<thisEnsScore ...
    
    smallData = mean(dataToUse(eligibleCells,:,:),3) - mean(mean(dataToUse(eligibleCells,bwin(1):bwin(2),:),2),3);
    smallCellVal = cellVal(eligibleCells);
    
    [s sidx] = sort(smallCellVal);
    
    figure(27);clf
    subplot(1,2,2);
    imagesc(smallData(fliplr(sidx),:))
    colorbar
    colormap rdbu
    % caxis([-0.75 0.75])
    caxis([-0.5 0.5])
    
    title(ens)
    
    subplot(1,2,1);
    r = rectangle('position',[6 -2 6 6]);
    r.FaceColor= [rgb('FireBrick') 0.5];
    r.LineStyle='none';
    
    hold on
    basVal = mean(mean(smallData(:,bwin(1):bwin(2)),2));
    fillPlot(smallData-basVal,[],'ci',rgb('steelblue'),'none',rgb('steelblue'),0.5);
    l = line([0 25],[0 0]);
    l.LineStyle=':';
    l.LineWidth=2;
    l.Color = rgb('grey');
    box off
    xlabel('Frames')
    ylabel('Mean Evoked dF/F')
    ylim([-0.1 0.1])
    
    pause
end

%% Plot just one
    %Nice ones 220302: 594 768 771 822 828 829 836 837 838 845 866 879 1005
    %1027 1030 1037 1053 1056 1058
ens = 838;%822 678; %578;% egL(i); %953;
thisEnsScore = outVars.popResponseEns(ens);


ind = outVars.ensIndNumber(ens)
h = outVars.ensHNumber(ens);
% sID = All(ind).out.exp.stimParams.Seq(h)
sID = h; %does this hold true if the first isn't a no hologram?

% [dfData, zdfData] =  computeDFFwithMovingBaseline(All(ind).out.exp.allData);
try 
    rawDat = All(ind).out.exp.dfData;
catch
    disp('err')
    if ~isfield(All(ind).out.exp,'dfData')
        disp('recomputing dfData')
        [dfData, zdfData] =  computeDFFwithMovingBaseline(All(ind).out.exp.allData);
        All(ind).out.exp.dfData = dfData;
        rawDat = dfData;
    end
end
% % % rawDat =dfData;% All(ind).out.exp.zdfData;
rawDatForCellVal = All(ind).out.exp.zdfData;

us = unique(All(ind).out.exp.stimID);

trialsToUse1 = All(ind).out.exp.lowMotionTrials &...
    All(ind).out.exp.lowRunTrials &...
    All(ind).out.exp.stimSuccessTrial &...
    All(ind).out.exp.stimID == us(sID) & ...
    All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);

dataToUse = rawDat(:,:,trialsToUse1); %order cells, frames, trials
win = All(ind).out.anal.recWinUsed;
bwin = All(ind).out.anal.bwinToUse;

cellVal = mean(mean(rawDatForCellVal(:,win(1):win(2),trialsToUse1),2),3)';
bVal = mean(mean(rawDatForCellVal(:,bwin(1):bwin(2),trialsToUse1),2),3)';
cellVal = cellVal-bVal;

trialsToUse2 = All(ind).out.exp.lowMotionTrials &...
    All(ind).out.exp.lowRunTrials &...
    All(ind).out.exp.stimSuccessTrial &...
    All(ind).out.exp.stimID == us(1) & ...
    All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);


noStimDataToUse = rawDat(:,:,trialsToUse2);


trialAverageData = mean(dataToUse,3);

offH = All(ind).out.exp.stimParams.roi{h};
thisOffTargetRisk = logical(All(ind).out.anal.offTargetRisk(offH,:));

eligibleCells  = ~thisOffTargetRisk...
    & ~All(ind).out.anal.ROIinArtifact';

%
% cellList = 1:numel(All(ind).out.anal.ROIinArtifact);
% eligibleCells  =  ~All(ind).out.anal.ROIinArtifact' ...
%     & ismember(cellList,All(ind).out.exp.holoTargets{offH});

smallData = mean(dataToUse(eligibleCells,:,:),3) - mean(mean(dataToUse(eligibleCells,bwin(1):bwin(2),:),2),3);
smallCellVal = cellVal(eligibleCells);

[s sidx] = sort(smallCellVal);

figure(27);clf
subplot(1,2,2);
imagesc(smallData(fliplr(sidx),:))
colorbar
colormap rdbu
% caxis([-0.75 0.75])
caxis([-0.6 0.6])

title(ens)

subplot(1,2,1);
r = rectangle('position',[6 -2 6 6]);
r.FaceColor= [rgb('FireBrick') 0.5];
r.LineStyle='none';

hold on
basVal = mean(mean(smallData(:,bwin(1):bwin(2)),2));
fillPlot(smallData-basVal,[],'ci',rgb('steelblue'),'none',rgb('steelblue'),0.5);
l = line([0 25],[0 0]);
l.LineStyle=':';
l.LineWidth=2;
l.Color = rgb('grey');
box off
xlabel('Frames')
ylabel('Mean Evoked dF/F')
ylim([-0.125 0.05])

%%
plotMask=1;

eligibleCellList = find(eligibleCells);
figure(28)
caxisLim =[-0.5 3];
colormap viridis
for i = (1:numel(smallCellVal)); %eligibleCellList)
    cellID = eligibleCellList(sidx(i)); %eligibleCellList(i);
    thisCellVal = smallCellVal(sidx(i));
    if 1;eligibleCells(cellID);
        clf;
        subplot(2,3,1)
        noStimDataToPlot = squeeze(noStimDataToUse(cellID,:,:))';
        imagesc(noStimDataToPlot)
        colorbar
        caxis(caxisLim)
        %     caxis([-3 3]);
        %     colormap rdbu
        title('No Stim')
        subplot(2,3,2)
        stimDatToPlot = squeeze(dataToUse(cellID,:,:))';
        imagesc(stimDatToPlot)
        colorbar
        caxis(caxisLim)
        
        FR = All(ind).out.info.FR; 
        xticks(0:FR:24)
        xticklabels(0:4);
        
        %     caxis([-3 3]);
        % %     colormap rdbu
        title(['Stim. ID: ' num2str(cellID)])
        
        subplot(2,3,3)
        fillPlot(noStimDataToPlot,[],'ci',rgb('DimGray'),'none',rgb('dimgray'),0.5);
        fillPlot(stimDatToPlot,[],'ci',rgb('Firebrick'),'none',rgb('firebrick'),0.5);
        title(['Val: ' num2str(thisCellVal)])
        
        subplot(2,3,6)
        %     fillPlot(noStimDataToPlot,[],'ci',rgb('DimGray'),'none',rgb('dimgray'),0.5);
        r = rectangle('position',[FR -2 6 FR]);
        r.FaceColor= [rgb('FireBrick') 0.5];
        r.LineStyle='none';
        
        hold on
        basVal = mean(mean(stimDatToPlot(:,bwin(1):bwin(2)),2));
        %         fillPlot(stimDatToPlot-basVal,[],'ci',rgb('steelblue'),'none',rgb('steelblue'),0.5);
        fillPlot(stimDatToPlot-basVal,[],'ci',rgb('darkmagenta'),'none',rgb('darkmagenta'),0.5);
%         fillPlot(stimDatToPlot-basVal,[],'ci',rgb('forestGreen'),'none',rgb('forestGreen'),0.5);
        
        l = line([0 25],[0 0]);
        l.LineStyle=':';
        l.LineWidth=2;
        l.Color = rgb('grey');
        box off
        xlabel('Time(s)')
        ylabel('Mean Evoked dF/F')
        ylim([-1.5 1])
                xticks(0:FR:24)
        xticklabels(0:4);
        
        if plotMask
            subplot(2,3,4)
            allDepth =All(ind).out.exp.allDepth;
            depth = allDepth(cellID);
            cellList = 1:numel(allDepth);
            %figure out the ID in as dat sees it
            IDinDat = numel(find(allDepth==depth & cellList'<=cellID));
            
            [dat, All(ind)] = pullDat(All(ind),depth);
            
            iscelldepth = [dat.stat.iscell];
            iscelldepth = find(iscelldepth);
            IDinDat = iscelldepth(IDinDat);
            
            sz = size(dat.mimg);
            blank= zeros(sz(1:2));
            
            xpix = dat.stat(IDinDat).xpix;
            ypix = dat.stat(IDinDat).ypix;
            lambda = dat.stat(IDinDat).lam;
            npix = dat.stat(IDinDat).npix;
            
            med = dat.stat(IDinDat).med;
            xySize = 7;
            
            for i=1:npix
                blank(xpix(i),ypix(i))=lambda(i);
            end
            sz = size(blank);
            
            blankSmall = blank(max(med(2)-xySize,1):min(med(2)+xySize,sz(1)),...
                max(med(1)-xySize,1):min(med(1)+xySize,sz(2)) );
            imagesc(blankSmall)
            axis square
            axis off
            colormap viridis
            
            muPerPx = 800/512;
            r = rectangle('Position',[1,xySize*2,10/muPerPx,1]);
            r.FaceColor= rgb('white');
            r.LineStyle = 'none';
        end
        %     depth = All(ind).out.exp.allDepth(cellID);
        %     [dat, All(ind)] = pullDat(All(ind),depth);
        %
        
        pause
    end
end

%% Plots in Figs
figure(2);clf
cell1 = 762;
cell2 = 706;

ens = 578;
ind = outVars.ensIndNumber(ens);

caxisLim =[-0.5 3];

% [dfData, zdfData] =  computeDFFwithMovingBaseline(All(ind).out.exp.allData);

rawDat =All(ind).out.exp.dfData;% All(ind).out.exp.zdfData;
rawDatForCellVal = All(ind).out.exp.zdfData;

us = unique(All(ind).out.exp.stimID);

trialsToUse1 = All(ind).out.exp.lowMotionTrials &...
    All(ind).out.exp.lowRunTrials &...
    All(ind).out.exp.stimSuccessTrial &...
    All(ind).out.exp.stimID == us(sID) & ...
    All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);

dataToUse = rawDat(:,:,trialsToUse1); %order cells, frames, trials

ax = subplot(2,4,1);
[All(ind)] = plotCellMask(All(ind),cell1,ax)

subplot(2,4,2)
stimDatToPlot = squeeze(dataToUse(cell1,:,:))';
imagesc(stimDatToPlot)
colorbar
caxis(caxisLim)

subplot(2,4,5);
r = rectangle('position',[6 -2 6 6]);
r.FaceColor= [rgb('FireBrick') 0.5];
r.LineStyle='none';

hold on
basVal = mean(mean(stimDatToPlot(:,bwin(1):bwin(2)),2));
fillPlot(stimDatToPlot-basVal,[],'ci',rgb('darkmagenta'),'none',rgb('darkmagenta'),0.5);

l = line([0 25],[0 0]);
l.LineStyle=':';
l.LineWidth=2;
l.Color = rgb('grey');
box off
xlabel('Frames')
ylabel('Mean Evoked dF/F')
ylim([-1 4])


ax = subplot(2,4,3);
[All(ind)] = plotCellMask(All(ind),cell2,ax)

subplot(2,4,4)
stimDatToPlot = squeeze(dataToUse(cell2,:,:))';
imagesc(stimDatToPlot)
colorbar
caxis(caxisLim)

subplot(2,4,7);
r = rectangle('position',[6 -2 6 6]);
r.FaceColor= [rgb('FireBrick') 0.5];
r.LineStyle='none';

hold on
basVal = mean(mean(stimDatToPlot(:,bwin(1):bwin(2)),2));
fillPlot(stimDatToPlot-basVal,[],'ci',rgb('steelblue'),'none',rgb('steelblue'),0.5);

l = line([0 25],[0 0]);
l.LineStyle=':';
l.LineWidth=2;
l.Color = rgb('grey');
box off
xlabel('Frames')
ylabel('Mean Evoked dF/F')
ylim([-2 2])