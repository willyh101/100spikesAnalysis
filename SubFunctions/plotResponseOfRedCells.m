function [outVars] = plotResponseOfRedCells(All,outVars,opts)
%% 
numExps = numel(All);

minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));
strtFrames = arrayfun(@(x) x.out.anal.recStartFrame,All);
totNumFrames = arrayfun(@(x) size(x.out.exp.zdfData,2),All);
minEndFrame = min(totNumFrames-minStrtFrame);

    clear popRespRed popRespNotRed popRespRedSEM popRespNotRedSEM
    clear popTSRed popTSNotRed nRedNotRed
c=0;
for ind = 1:numExps
    
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    
    us = unique(All(ind).out.exp.stimID);
    
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;
    
    isRed = All(ind).out.red.isRed;
    
    v=1;
        
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    
    for i= 1:numStims
        holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
        c=c+1;
        if i==1
            cellsToUse = ~ROIinArtifact';
        else
            cellsToUse = ~ROIinArtifact'  & ~offTargetRisk(holo,:);
        end
        popRespRed(c) = nanmean(squeeze(respMat(i,v,cellsToUse & isRed) -...
            baseMat(i,v,cellsToUse & isRed)));
        popRespRedSEM(c) = sem(squeeze(respMat(i,v,cellsToUse & isRed) -...
            baseMat(i,v,cellsToUse & isRed)));
        popRespNotRed(c) = mean(squeeze(respMat(i,v,cellsToUse & ~isRed) -...
            baseMat(i,v,cellsToUse & ~isRed)));
        popRespNotRedSEM(c) = sem(squeeze(respMat(i,v,cellsToUse & ~isRed) -...
            baseMat(i,v,cellsToUse & ~isRed)));
        isControl(c) = i==1;
        
        %TS code
        strtFrame = All(ind).out.anal.recStartFrame;
        newStart = strtFrame-minStrtFrame+1;
        newEnd = newStart+minEndFrame-1;
        s = us(i);
        
        
        dat = All(ind).out.exp.zdfData(cellsToUse & isRed,newStart:newEnd,trialsToUse &...
            All(ind).out.exp.stimID==s &...
            All(ind).out.exp.visID==v );
        popTSRed(c,:) = mean(mean(dat,3),1);
        dat = All(ind).out.exp.zdfData(cellsToUse & ~isRed,newStart:newEnd,trialsToUse &...
            All(ind).out.exp.stimID==s &...
            All(ind).out.exp.visID==v );
        popTSNotRed(c,:) = mean(mean(dat,3),1);
        nRedNotRed(c,1) = sum(cellsToUse & isRed);
        nRedNotRed(c,2) = sum(cellsToUse & ~isRed);
    end
 

 
end

outVars.popRespRed          = popRespRed;
outVars.popRespRedSEM       = popRespRedSEM;
outVars.popRespNotRed       = popRespNotRed;
outVars.popRespNotRedSEM    = popRespNotRedSEM;
outVars.popTSRed            = popTSRed;
outVars.popTSNotRed         = popTSNotRed;

 ensToPlot =outVars.ensemblesToUse;


%  figure(23);clf
%  plot(popRespRed,'o')
%  hold on
%  plot(popRespNotRed,'o')
% 
%  figure(23);clf
%  plot(popRespRed,popRespNotRed,'o')
%  refline(1)
 
% Means Section
popRespRed(isControl)=[];
popRespNotRed(isControl)=[];


figure(233);clf
% sp3 = scatter(ones(1,numel(popRespRed(ensToPlot))), popRespRed(ensToPlot), 'filled');
% sp3.MarkerEdgeColor = 'b';%'r'
% sp3.MarkerFaceAlpha = 0.6;
% sp3.MarkerEdgeColor = 'none';
% hold on
% e3 = errorbar(1, mean(popRespRed(ensToPlot)), sem(popRespRed(ensToPlot)));
% e3.LineWidth = 2;
% e3.Color = 'k';
% 
% sp32 = scatter(ones(1,numel(popRespNotRed(ensToPlot)))+1, popRespNotRed(ensToPlot), 'filled');
% sp3.MarkerEdgeColor = 'b';%'r'
% sp3.MarkerFaceAlpha = 0.6;
% sp3.MarkerEdgeColor = 'none';
% e32 = errorbar(2, mean(popRespNotRed(ensToPlot)), sem(popRespNotRed(ensToPlot)));
% e32.LineWidth = 2;
% e32.Color = 'k';
% xlim([0.5 2.5])
p = fancyPlotSpread([popRespRed(ensToPlot); popRespNotRed(ensToPlot)]', {opts.redCellName, 'Pyramidal Cells'});
ylabel('Mean Population Response')
% xlim([0.7 2.3])


figure(23);clf
sp = scatter(outVars.numCellsEachEns(ensToPlot),popRespRed(ensToPlot),'filled');
sp.MarkerFaceAlpha = 0.6;
sp.MarkerEdgeColor = 'none';
hold on
sp2=scatter(outVars.numCellsEachEns(ensToPlot)+0.75,popRespNotRed(ensToPlot),'filled');
sp2.MarkerFaceAlpha = 0.6;
sp2.MarkerEdgeColor = 'none';
sp2.MarkerFaceColor = rgb('tomato');


  
  
 r = refline(0);
 r.LineStyle=':';
 r.Color=rgb('grey');
 r.LineWidth=2;
 ylabel('Population Response')
 xlabel('Ensemble Size')

 %now determine statistical signifigance
 ensCat = unique(outVars.numCellsEachEns(ensToPlot));
 clear redPVal redGrandMean redGrandSEM notRedGrandMean notRedGrandSEM
 for i =1:numel(ensCat);
     redPVal(i) = signrank(popRespRed(ensToPlot & outVars.numCellsEachEns==ensCat(i)),...
         popRespNotRed(ensToPlot & outVars.numCellsEachEns==ensCat(i)) );
     redGrandMean(i) = nanmean(popRespRed(ensToPlot & outVars.numCellsEachEns==ensCat(i)));
     redGrandSEM(i) = sem(popRespRed(ensToPlot & outVars.numCellsEachEns==ensCat(i)));
     
     notRedGrandMean(i) = nanmean(popRespNotRed(ensToPlot & outVars.numCellsEachEns==ensCat(i)));
     notRedGrandSEM(i) = sem(popRespNotRed(ensToPlot & outVars.numCellsEachEns==ensCat(i)));
     
     disp(['Size ' num2str(ensCat(i)) ' Red mean: ' num2str(redGrandMean(i)) '. Not Red: ' num2str(notRedGrandMean(i)) ]);
     disp(['Singed Rank: ' num2str(redPVal(i))])
 end
 e = errorbar(ensCat+0.25,redGrandMean,redGrandSEM);
 e.LineWidth = 2;
 e.Color =rgb('royalblue'); %'r'
 e = errorbar(ensCat+0.5,notRedGrandMean,notRedGrandSEM);
 e.LineWidth = 2;
 e.Color =rgb('tomato');
 
%  legend([sp(1) sp2(1)],'Red Cells', 'Not Red Cells')
 legend([sp(1) sp2(1)], opts.redCellName, 'Pyramidal Cells')

 %TS Section
 popTSRed(isControl,:)=[];
 popTSNotRed(isControl,:)=[];
 
 baseline = 1;
 if baseline
     popTSRed = popTSRed - mean(popTSRed(:,1:minStrtFrame),2);
     popTSNotRed = popTSNotRed - mean(popTSNotRed(:,1:minStrtFrame),2);
 end
 
 colorLim = [-0.5 0.5];
 figure(22);clf
 ax1(1) = subplot(2,2,1);
 imagesc(popTSRed(ensToPlot,:))
 colormap rdbu
 caxis(colorLim)
 title('Red Cells')
 ylabel('Ensemble')
 xlabel('Frame')
 
 ax(1) = subplot(2,2,3);
 fillPlot(popTSRed(ensToPlot,:),[],'ci');
 
 ax1(2) = subplot(2,2,2);
  imagesc(popTSNotRed(ensToPlot,:))
 colormap rdbu
 caxis(colorLim)
 title('Not Red Cells')
 xlabel('Frame')
 
 ax(2) = subplot(2,2,4);
  fillPlot(popTSNotRed(ensToPlot,:),[],'ci');
linkaxes(ax);

figure(24);
clf;

rc =rectangle('position',[minStrtFrame -0.04 6 0.07]);
rc.FaceColor = [rgb('FireBrick') 0.25];
rc.LineStyle = 'none';

hold on
lineCol = rgb('SteelBlue');% rgb('FireBrick');
edgeCol = 'none';
faceCol = rgb('SteelBlue');% rgb('FireBrick');
faceAlpha = 0.5;
fp1 =  fillPlot(popTSRed(ensToPlot,:),[],'ci',lineCol,edgeCol,faceCol,faceAlpha);
lineCol = rgb('Black');
edgeCol = 'none';
faceCol = rgb('DimGray');
faceAlpha = 0.5;
 fp2 = fillPlot(popTSNotRed(ensToPlot,:),[],'ci',lineCol,edgeCol,faceCol,faceAlpha);
 legend([fp1(1) fp2(1)], opts.redCellName,'Pyramidal Cells')
 
 xlabel('Frame')
 ylabel('\DeltaZ-Score dF/F')
 
 