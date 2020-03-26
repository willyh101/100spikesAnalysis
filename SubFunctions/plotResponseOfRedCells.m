function [outVars] = plotResponseOfRedCells(All,outVars,opts)
%% 
numExps = numel(All);
minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));


    clear popRespRed popRespNotRed popRespRedSEM popRespNotRedSEM
    clear popTSRed popTSNotRed nRedNotRed
c=0;
for ind = 1:numExps;
    
    
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
            if i==1;
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
            s = us(i); 
            
            
           dat = All(ind).out.exp.zdfData(cellsToUse & isRed,newStart:end,trialsToUse &...
                All(ind).out.exp.stimID==s &...
                All(ind).out.exp.visID==v );
            popTSRed(c,:) = mean(mean(dat,3),1);
            dat = All(ind).out.exp.zdfData(cellsToUse & ~isRed,newStart:end,trialsToUse &...
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


 figure(23);clf
 sp = scatter(outVars.numCellsEachEns(ensToPlot),popRespRed(ensToPlot),'o');
 sp.MarkerEdgeColor = 'r';
 hold on
 sp=scatter(outVars.numCellsEachEns(ensToPlot)+0.5,popRespNotRed(ensToPlot),'o');
  sp.MarkerEdgeColor = 'k';

 r = refline(0);
 r.LineStyle=':';
 r.Color=rgb('grey');
 r.LineWidth=2;
 ylabel('Population Response')
 xlabel('Ensemble Size')
 legend({'Red Cells'; 'Not Red Cells'})

 
 %TS Section
 popTSRed(isControl,:)=[];
 popTSNotRed(isControl,:)=[];
 
 baseline = 1;
 if baseline
     popTSRed = popTSRed - mean(popTSRed(:,1:minStrtFrame),2);
     popTSNotRed = popTSNotRed - mean(popTSNotRed(:,1:minStrtFrame),2);
 end
 figure(22);clf
 subplot(1,2,1)
 imagesc(popTSRed(ensToPlot,:))
 colormap rdbu
 caxis([-0.25 0.25])
 title('Red Cells')
 ylabel('Ensemble')
 xlabel('Frame')
 
 subplot(1,2,2)
  imagesc(popTSNotRed(ensToPlot,:))
 colormap rdbu
 caxis([-0.25 0.25])
 title('Not Red Cells')
 xlabel('Frame')
 
 
 