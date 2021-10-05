%% stand alone plot of some example Data

meanEnsResp = mean(outVars.popResponseEns(ensemblesToUse));
stdEnsResp = std(outVars.popResponseEns(ensemblesToUse));

eligibleEns = outVars.popResponseEns<meanEnsResp & outVars.popResponseEns<meanEnsResp>(meanEnsResp-stdEnsResp) & ensemblesToUse';


ens = 953;
thisEnsScore = outVars.popResponseEns(ens);


ind = outVars.ensIndNumber(ens);
h = outVars.ensHNumber(ens); 
sID = h+1; %does this hold true if the first isn't a no hologram?


rawDat = All(ind).out.exp.zdfData;
rawDatForCellVal = All(ind).out.exp.zdfData;

us = All(ind).out.exp.uniqueStims;

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


thisOffTargetRisk = logical(All(ind).out.anal.offTargetRisk(h,:));

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
%%
eligibleCellList = find(eligibleCells);
figure(28)
colormap viridis
for i = 1:numel(eligibleCellList)
    cellID = eligibleCellList(i);
    clf;
    subplot(2,3,1)
    noStimDataToPlot = squeeze(noStimDataToUse(eligibleCellList(i),:,:))';
    imagesc(noStimDataToPlot)
    colorbar
    caxis([0.5 3])
%     caxis([-3 3]);
%     colormap rdbu
    title('No Stim')
    subplot(2,3,2)
    stimDatToPlot = squeeze(dataToUse(eligibleCellList(i),:,:))';
    imagesc(stimDatToPlot)
    colorbar
        caxis([0.5 3])

%     caxis([-3 3]);
% %     colormap rdbu
    title(['Stim. ID: ' num2str(cellID)])
    
    subplot(2,3,3)
    fillPlot(noStimDataToPlot,[],'ci',rgb('DimGray'),'none',rgb('dimgray'),0.5);
    fillPlot(stimDatToPlot,[],'ci',rgb('Firebrick'),'none',rgb('firebrick'),0.5);
    pause
end
    