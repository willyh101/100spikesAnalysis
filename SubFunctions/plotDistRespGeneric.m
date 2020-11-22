function [eHandle] = plotDistRespGeneric(RespToPlot,outVars,opts,axesHandle);

% popResponseDist = outVars.popResponseDist;
% numSpikesEachStim = outVars.numSpikesEachStim;
% popResponseNumCells =outVars.popResponseNumCells;
% 
ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;

distBins = opts.distBins;

xaxisrange = opts.distAxisRange;

%% Plot Pop Response by Distance
% popDistAll = cell2mat(RespToPlot');
% popDistAll = cell2mat(popResponseDistVis');

popDist =RespToPlot;% popDistAll(numSpikesEachStim~=0,:);
% 
% popNumCells = cell2mat(popResponseNumCells');
% popNCells = popNumCells(numSpikesEachStim~=0,:);


ensSizes = unique(numCellsEachEns(ensemblesToUse))   ;

numEnsembles = numel(ensSizes);
if numEnsembles == 1
    colorList{1} = rgb('black');
else
colorList = colorMapPicker(numEnsembles,outVars.defaultColorMap);
end
% if numEnsembles==3
%     colorList = {rgb('DarkBlue') rgb('steelblue') rgb('gold')};
% else
%     cl=[];
%     colorList = [];%{rgb('gray')};
%     cl = colormap('plasma');%colormap('rdbu');
% %     nColors = 3; numel(MRDat);
%     ncs = round(linspace(1,size(cl,1),numEnsembles));
%     for i=1:numEnsembles;
%         colorList{end+1} =  cl(ncs(i),:);
%     end
% end
% figure(9);clf
subplot(axesHandle);
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat = popDist(ensemblesToUse & numCellsEachEns==ensSizes(i),:);
meanDat = nanmean(dat);
stdDat = nanstd(dat);
numpDat = sum(~isnan(dat));
semDat = stdDat./sqrt(numpDat);


hold on
distBinSize = distBins(2)-distBins(1);
eHandle{i} = errorbar(distBins(2:end)-distBinSize/2,meanDat,semDat,'linewidth',2,'color',colorList{i});
end
r = refline(0);
r.LineStyle=':';
r.Color = rgb('grey');
r.LineWidth = 2;
xlabel('Distance from a target')
ylabel('Population Response (mean of ensembles'' pop response)')
xlim(xaxisrange)
if numEnsembles ==3;
    legend('Small', 'Medium', 'Big')
else
    legend(cellfun(@num2str,num2cell(ensSizes),'uniformoutput',0))
end