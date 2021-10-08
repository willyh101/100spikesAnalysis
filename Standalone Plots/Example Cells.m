%% stand alone plot of some example Data

meanEnsResp = mean(outVars.popResponseEns(ensemblesToUse));
stdEnsResp = std(outVars.popResponseEns(ensemblesToUse));

eligibleEns = outVars.popResponseEns<meanEnsResp & outVars.popResponseEns<meanEnsResp>(meanEnsResp-stdEnsResp) & ensemblesToUse';

egL = find(eligibleEns);
for i =1:sum(eligibleEns)
ens = egL(i); %953;
thisEnsScore = outVars.popResponseEns(ens);


ind = outVars.ensIndNumber(ens);
h = outVars.ensHNumber(ens);
% sID = All(ind).out.exp.stimParams.Seq(h)
sID = h; %does this hold true if the first isn't a no hologram?


rawDat = All(ind).out.exp.zdfData;
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
imagesc(smallData(fliplr(sidx),:))
colorbar
colormap rdbu
caxis([-0.75 0.75])
title(ens)
pause
end

%% Plot just one

ens = 481; egL(i); %953;
thisEnsScore = outVars.popResponseEns(ens);


ind = outVars.ensIndNumber(ens);
h = outVars.ensHNumber(ens);
% sID = All(ind).out.exp.stimParams.Seq(h)
sID = h; %does this hold true if the first isn't a no hologram?


rawDat = All(ind).out.exp.zdfData;
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
imagesc(smallData(fliplr(sidx),:))
colorbar
colormap rdbu
caxis([-0.75 0.75])
title(ens)
%%
eligibleCellList = find(eligibleCells);
figure(28)
caxisLim =[-0.5 3];
colormap viridis
for i = 1:numel(smallCellVal); %eligibleCellList)
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

%     caxis([-3 3]);
% %     colormap rdbu
    title(['Stim. ID: ' num2str(cellID)])
    
    subplot(2,3,3)
    fillPlot(noStimDataToPlot,[],'ci',rgb('DimGray'),'none',rgb('dimgray'),0.5);
    fillPlot(stimDatToPlot,[],'ci',rgb('Firebrick'),'none',rgb('firebrick'),0.5);
    title(['Val: ' num2str(thisCellVal)])
    
        subplot(2,3,6)
%     fillPlot(noStimDataToPlot,[],'ci',rgb('DimGray'),'none',rgb('dimgray'),0.5);
    fillPlot(stimDatToPlot,[],'ci',rgb('Firebrick'),'none',rgb('firebrick'),0.5);
    l = line([0 25],[0 0]);
    l.LineStyle=':';
    pause
    end
end
    