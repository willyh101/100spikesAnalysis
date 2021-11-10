function plotCellMean(All,outVars,ens,cell,color,ax)

ind = outVars.ensIndNumber(ens);
sID = outVars.ensHNumber(ens);

rawDat =All(ind).out.exp.dataToUse;% All(ind).out.exp.zdfData;

us = unique(All(ind).out.exp.stimID);

trialsToUse1 = All(ind).out.exp.lowMotionTrials &...
    All(ind).out.exp.lowRunTrials &...
    All(ind).out.exp.stimSuccessTrial &...
    All(ind).out.exp.stimID == us(sID) & ...
    All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);

dataToUse = rawDat(:,:,trialsToUse1); %order cells, frames, trials

stimDatToPlot = squeeze(dataToUse(cell,:,:))';

subplot(ax)
r = rectangle('position',[6 -2 6 6]);
r.FaceColor= [rgb('FireBrick') 0.5];
r.LineStyle='none';

win = All(ind).out.anal.recWinUsed;
bwin = All(ind).out.anal.bwinToUse;

hold on
basVal = mean(mean(stimDatToPlot(:,bwin(1):bwin(2)),2));
fillPlot(stimDatToPlot-basVal,[],'ci',color,'none',color,0.5);

l = line([0 25],[0 0]);
l.LineStyle=':';
l.LineWidth=2;
l.Color = rgb('grey');
box off
xlabel('Frames')
ylabel('Mean Evoked dF/F')
ylim([-1 4])

FR = All(ind).out.info.FR;

ticVals = [1:FR:size(rawDat,2)];
xticks(ticVals)
xticklabels(0:numel(ticVals))