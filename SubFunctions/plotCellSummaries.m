function plotCellSummaries(All,outVars,ens,cell1,figNum)

figure(figNum);clf

ind = outVars.ensIndNumber(ens);
sID = outVars.ensHNumber(ens);
caxisLim = [-0.25 3];

rawDat =All(ind).out.exp.dfData;% All(ind).out.exp.zdfData;
rawDatForCellVal = All(ind).out.exp.zdfData;

us = unique(All(ind).out.exp.stimID);

trialsToUse1 = All(ind).out.exp.lowMotionTrials &...
    All(ind).out.exp.lowRunTrials &...
    All(ind).out.exp.stimSuccessTrial &...
    All(ind).out.exp.stimID == us(sID) & ...
    All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);

dataToUse = rawDat(:,:,trialsToUse1); %order cells, frames, trials

ax = subplot(2,2,1);
[All(ind)] = plotCellMask(All(ind),cell1,ax);

subplot(2,2,2)
stimDatToPlot = squeeze(dataToUse(cell1,:,:))';
imagesc(stimDatToPlot)
colorbar
caxis(caxisLim)

subplot(2,2,3);
r = rectangle('position',[6 -2 6 6]);
r.FaceColor= [rgb('FireBrick') 0.5];
r.LineStyle='none';

win = All(ind).out.anal.recWinUsed;
bwin = All(ind).out.anal.bwinToUse;

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
