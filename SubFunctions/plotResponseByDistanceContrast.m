function plotResponseByDistanceContrast(outVars,opts)
popResponseAllDist = outVars.popResponseAllDist;
popResponseAllDistSub = outVars. popResponseAllDistSub;
popResponseAllDistSubVis = outVars.popResponseAllDistSubVis;

numSpikesEachStim = outVars.numSpikesEachStim;
popResponseNumCells =outVars.popResponseNumCells;

ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;

distBins = opts.distBins;
%% as above but for full and no vis Conditions

figure(11);clf


    
popDatNoVis = cell2mat(cellfun(@(x) permute(x(2:end,1,:),[1 3 2]), popResponseAllDist,'uniformoutput',0)');
% popDatNoVisNoStim = squeeze(cell2mat(cellfun(@(x) (x(1,1,:)), popResponseAllDist,'uniformoutput',0)'));

popDatMaxVis = cell2mat(cellfun(@(x) permute(x(2:end,end,:),[1 3 2]), popResponseAllDist,'uniformoutput',0)');
% popDatMaxVisNoStim = squeeze(cell2mat(cellfun(@(x) (x(1,end,:)), popResponseAllDist,'uniformoutput',0)'));
popDatMaxVisSubtracted = cell2mat(cellfun(@(x) permute(x(2:end,end,:),[1 3 2]), popResponseAllDistSub,'uniformoutput',0)');
popDatMaxVisSubVis = cell2mat(cellfun(@(x) permute(x(2:end,end,:),[1 3 2]), popResponseAllDistSubVis,'uniformoutput',0)');


divider = 1;
popDatVis2 = cell2mat(cellfun(@(x) squeeze(x(2:end,round(size(x,2)/divider),:)), popResponseAllDist,'uniformoutput',0)');
popDatVis2Subtracted = cell2mat(cellfun(@(x) squeeze(x(2:end,round(size(x,2)/divider),:)), popResponseAllDistSub,'uniformoutput',0)');


ensSizes = unique(numCellsEachEns(ensemblesToUse))   ;
numEnsembles = numel(ensSizes);

if numEnsembles==3
    colorList = {rgb('DarkBlue') rgb('steelblue') rgb('gold')};
else
    cl=[];
    colorList = [];%{rgb('gray')};
    cl = colormap('plasma');%colormap('rdbu');
%     nColors = 3; numel(MRDat);
    ncs = round(linspace(1,size(cl,1),numEnsembles));
    for i=1:numEnsembles;
        colorList{end+1} =  cl(ncs(i),:);
    end
end



subplot(1,2,1)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatNoVis(ensemblesToUse & numCellsEachEns==ensSizes(i) ,:);
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
xlim([0 350])
if numEnsembles ==3;
    legend('Small', 'Medium', 'Big')
else
    legend(cellfun(@num2str,num2cell(ensSizes),'uniformoutput',0))
end
title('Holographic induced changes 0 Contrast')

subplot(1,2,2)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatMaxVisSubtracted(ensemblesToUse & numCellsEachEns==ensSizes(i) ,:);
% dat = popDatMaxVis(ensemblesToUse & numCellsEachEns==ensSizes(i) & highVisPercentInd ,:);

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
xlim([0 350])
if numEnsembles ==3;
    legend('Small', 'Medium', 'Big')
else
    legend(cellfun(@num2str,num2cell(ensSizes),'uniformoutput',0))
end
title('Holographic induced changes Max Contrast')


figure(12);clf
contrastsToView = [6 3 2 1.5 1.25 1] ;%I know its weird i just wanted to be able to catch times that we were using different contrasts, will work out to 1:6 if there are 6 contrasts; 1:6;
for c=1:numel(contrastsToView)
ax(c) = subplot(1,numel(contrastsToView),c);
divider = contrastsToView(c);
popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistSub,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDist,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistSubVis,'uniformoutput',0)');
% popDatToPlot = cell2mat(cellfun(@(x) squeeze(x(2:end,max(round(size(x,2)/divider),1),:)), popResponseAllDistVis,'uniformoutput',0)');

for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDatToPlot(ensemblesToUse & numCellsEachEns==ensSizes(i) ,:);
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
xlim([0 400])
if numEnsembles ==3;
    legend('Small', 'Medium', 'Big')
else
    legend(cellfun(@num2str,num2cell(ensSizes),'uniformoutput',0))
end
title(['Contrast ' num2str(c) ] );



end

linkaxes([ax(:)])