function [All outVars] = createTSPlot(All,outVars)

%% Create time series plot

numExps = numel(All);
ensemblesToUse = outVars.ensemblesToUse; 
IndsUsed = outVars.IndsUsed;

minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));

clear allMeanTS allStdTS allnumTS allMeanTSVis allStdTSVis allnumTSVis
for ind=1:numExps
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    strtFrame = All(ind).out.anal.recStartFrame;
    newStart = strtFrame-minStrtFrame+1;
    
    pVisR = All(ind).out.anal.pVisR;
    
    clear mRespTS sRespTS nResp
    for i = 1:numel(us)
        s = us(i);
        h = All(ind).out.exp.stimParams.roi{i};
        
        if h>0
            tg = All(ind).out.exp.holoTargets{h};
            tg(isnan(tg))=[];
        else
            tg=[];
        end
        cellList = 1:numel(ROIinArtifact);
        
        if i==1
            cellsToUse = ~ROIinArtifact';% & pVisR<0.05 ;
        else
            cellsToUse = ~ROIinArtifact' &...
                ~offTargetRisk(h,:) &...
                ~ismember(cellList,tg);% & pVisR<0.05;
        end
        
        for k=1:numel(vs)
            v=vs(k);
            
            dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end,trialsToUse &...
                All(ind).out.exp.stimID==s &...
                All(ind).out.exp.visID==v );
            
            mDat = mean(dat,3);
            mmDat = mean(mDat,1); %pop Average
            sdDat = std(mDat);
            nDat = size(mDat,1);
            
            
            mRespTS(i,:,k) = mmDat; % mean response time series
            sRespTS(i,:,k) = sdDat; % std response time series (by cell);
            nResp(i,k) = nDat;
        end
        
    end
    All(ind).out.anal.mRespTS= mRespTS;
    All(ind).out.anal.sRespTS= sRespTS;
    All(ind).out.anal.nResp = nResp;
    
    
    allMeanTS{ind} = mRespTS(:,:,1);
    allStdTS{ind} = sRespTS(:,:,1);
    allnumTS{ind} = nResp(:,1);
    
    allMeanTSVis{ind} = mRespTS;
    allStdTSVis{ind} = sRespTS;
    allnumTSVis{ind} = nResp;
    
end

baseline=1;

shortestRec = min(cellfun(@(x) size(x,2),allMeanTS));
temp = cellfun(@(x) x(2:end,1:shortestRec),allMeanTS,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
meanTSSquare = cat(1,temp{:});

%no Resp
temp = cellfun(@(x) x(1,1:shortestRec),allMeanTS,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
meanTSSquareNR = cat(1,temp{:});

temp = cellfun(@(x) x(2:end,1:shortestRec),allStdTS,'uniformoutput',0);
stdTSSquare = cat(1,temp{:});

figure(25);clf
subplot(2,2,2)
imagesc(meanTSSquare(ensemblesToUse,:))
% imagesc(meanTSSquare(:,:))
title('stim')

colormap rdbu
caxis([-0.2 0.2])

clear ax
ax(1)=subplot(2,2,4);
fillPlot(meanTSSquare(ensemblesToUse,:),[],'ci');

subplot(2,2,1)
imagesc(meanTSSquareNR(IndsUsed,:))
title('NoStim')

colormap rdbu
caxis([-0.1 0.1])
ax(2) = subplot(2,2,3);
if numel(IndsUsed)>1
fillPlot(meanTSSquareNR(IndsUsed,:),[],'ci');
else
    p= plot(meanTSSquareNR(IndsUsed,:));
    p.Color = 'k';
    p.LineWidth = 2;
end
linkaxes(ax);