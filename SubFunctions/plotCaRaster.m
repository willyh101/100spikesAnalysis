function plotCaRaster(All,outVars,ens,cell,ax)
caxisLim = [-0.25 3];

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


subplot(ax)
stimDatToPlot = squeeze(dataToUse(cell,:,:))';
imagesc(stimDatToPlot)
colorbar
caxis(caxisLim)

FR = All(ind).out.info.FR;

ticVals = [1:FR:size(rawDat,2)];
xticks(ticVals)
xticklabels(0:numel(ticVals))