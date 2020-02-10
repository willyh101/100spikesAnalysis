%% Ian Iso Tuned Trials
degrees = [nan 0    45    90   135   180   225   270   315];
for ind=1:numExps
    EnsPref = All(ind).out.anal.ensemblePrefDeg;    
    us = unique(All(ind).out.exp.stimID);
    stimID = All(ind).out.exp.stimID;
    try
        visDeg = All(ind).out.exp.visCond(2,:);
    catch
        visDeg = nan(size(stimID));
    end
    
    ortho1 = EnsPref-90;
    ortho1(ortho1<0) = ortho1(ortho1<0)+360;
    
    ortho2 = EnsPref+90;
    ortho2(ortho2>360) = ortho2(ortho2>360)-360;
    
    clear isoTunedTrial orthoTunedTrial
    for i=1:numel(stimID)
        s=stimID(i);
        sidx = find(us==s);
        holo = All(ind).out.exp.stimParams.roi{sidx};
        
        thisVisDeg = visDeg(i);
        if holo>0
            thisStimDeg = EnsPref(holo);
            thisStimOrtho1 =ortho1(holo);
            thisStimOrtho2 = ortho2(holo);
        else
            thisStimDeg=nan;
            thisStimOrtho2 = nan;
            thisStimOrtho1 = nan;
        end
         
        isoTunedTrial(i) = thisVisDeg==thisStimDeg;
        orthoTunedTrial(i) = thisVisDeg==thisStimOrtho1 | thisVisDeg==thisStimOrtho2;
    end
    All(ind).out.anal.isoTunedTrial = isoTunedTrial;
    All(ind).out.anal.orthoTunedTrial = orthoTunedTrial;
    countIsoTrial(ind) = sum(isoTunedTrial);
    countOrthoTrial(ind) = sum(orthoTunedTrial);
end
        
%% Create time series plot
minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));

visAlphaAnalysis = 0.05;

clear allMeanTS allStdTS allnumTS allMeanTSVis allStdTSVis allnumTSVis
for ind=1:numExps
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial ;
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    strtFrame = All(ind).out.anal.recStartFrame;
    newStart = strtFrame-minStrtFrame+1;
    
    pVisR = All(ind).out.anal.pVisR;
    
    clear mRespTS sRespTS nResp
    mRespIso=[]; mRespIsoSub=[];
    sRespIso=[];
    nRespIso=[];
    mRespOrtho=[]; mRespOrthoSub=[];
    sRespOrtho=[];
    nRespOrtho=[];
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
            cellsToUse = ~ROIinArtifact' & pVisR<visAlphaAnalysis ;
        else
            thisPref = All(ind).out.anal.ensemblePrefDeg(h);
            cellOri = All(ind).out.anal.prefOri;
            cellOri = degrees(cellOri);
            
            cellsToUse = ~ROIinArtifact'  &...
                ~offTargetRisk(h,:) &...
                ~ismember(cellList,tg) &...
                pVisR<visAlphaAnalysis &... 
                cellOri == thisPref;
        end
        
        for k=1:numel(vs)
            v=vs(k);
            
            thisTrialsToUse = trialsToUse &...
                All(ind).out.exp.stimID==s &...
                All(ind).out.exp.visID==v &...
                All(ind).out.anal.isoTunedTrial;
 
            dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end, thisTrialsToUse);
            
            mDat = mean(dat,3);
            mmDat = mean(mDat,1); %pop Average
            sdDat = std(mDat);
            nDat = size(mDat,1);
            mRespIso(i,:,k) = mmDat; % mean response time series
            sRespIso(i,:,k) = sdDat; % std response time series (by cell);
            nRespIso(i,k) = nDat;
            
             thisTrialsToUse = trialsToUse &...
                All(ind).out.exp.stimID==us(1) &...
                All(ind).out.exp.visID==v;
            dat2 = All(ind).out.exp.zdfData(cellsToUse,newStart:end, thisTrialsToUse);
            mDat2 = mean(dat2,3);
            mmDat2 = mean(mDat2,1); %pop Average
            sdDat2 = std(mDat2);
            nDat2 = size(mDat2,1);
            
            mRespIsoSub(i,:,k) = mmDat-mmDat2; % mean response time series

            
            
%             if sum(thisTrialsToUse)>0
%                 mRespIso(i,:) = mmDat;
%                 sRespIso(i,:) = sdDat;
%                 nRespIso(i,:) = nDat;
%             end
            
            thisTrialsToUse = trialsToUse &...
                All(ind).out.exp.stimID==s &...
                All(ind).out.exp.visID==v &...
                All(ind).out.anal.orthoTunedTrial;
             dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end, thisTrialsToUse);
            
            mDat = mean(dat,3);
            mmDat = mean(mDat,1); %pop Average
            sdDat = std(mDat);
            nDat = size(mDat,1);
            
            mRespOrtho(i,:,k) = mmDat;
            sRespOrtho(i,:,k) = sdDat;
            nRespOrtho(i,:,k) = nDat;
            mRespOrthoSub(i,:,k) = mmDat-mmDat2; % mean response time series

            
           
        end
        
    end
    
%     All(ind).out.anal.mRespTS= mRespTS;
%     All(ind).out.anal.sRespTS= sRespTS;
%     All(ind).out.anal.nResp = nResp;
%     
%     
%     allMeanTS{ind} = mRespTS(:,:,1);
%     allStdTS{ind} = sRespTS(:,:,1);
%     allnumTS{ind} = nResp(:,1);
%     
%     allMeanTSVis{ind} = mRespTS;
%     allStdTSVis{ind} = sRespTS;
%     allnumTSVis{ind} = nResp;
    
allmRespIso{ind}=nanmean(mRespIso,3);
allsRespIso{ind}=nanmean(sRespIso,3);
allnRespIso{ind}=nanmean(nRespIso,2);

allmRespOrtho{ind}=nanmean(mRespOrtho,3);
allsRespOrtho{ind}=nanmean(sRespOrtho,3);
allnRespOrtho{ind}=nanmean(nRespOrtho,2);

allmRespIsoSub{ind}=nanmean(mRespIsoSub,3);
allmRespOrthoSub{ind}=nanmean(mRespOrthoSub,3);


end
%%
baseline=1;

shortestRec = min(cellfun(@(x) size(x,2),allmRespIso));

temp = cellfun(@(x) x(2:end,1:shortestRec),allmRespIso,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
allmRespIsoSquare = cat(1,temp{:});
sum(~isnan(allmRespIsoSquare(ensemblesToUse,1)))


temp = cellfun(@(x) x(2:end,1:shortestRec),allmRespOrtho,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
allmRespOrthoSquare = cat(1,temp{:});
sum(~isnan(allmRespOrthoSquare(ensemblesToUse,1)))

%%
baseline=1;

shortestRec = min(cellfun(@(x) size(x,2),allmRespIsoSub));

temp = cellfun(@(x) x(2:end,1:shortestRec),allmRespIsoSub,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
allmRespIsoSquare = cat(1,temp{:});
sum(~isnan(allmRespIsoSquare(ensemblesToUse,1)))


temp = cellfun(@(x) x(2:end,1:shortestRec),allmRespOrthoSub,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
allmRespOrthoSquare = cat(1,temp{:});
sum(~isnan(allmRespOrthoSquare(ensemblesToUse,1)))

%%
figure(25);clf
subplot(2,3,2)
 valid = ~isnan(allmRespIsoSquare(:,1))' & ensemblesToUse;
 allmRespIsoSquare(~valid,:)=[];

 plot(allmRespIsoSquare')
 subplot(2,3,3)
 valid = ~isnan(allmRespOrthoSquare(:,1))' & ensemblesToUse;
 allmRespOrthoSquare(~valid,:)=[];

 plot(allmRespOrthoSquare')

%%
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
imagesc(allmRespIsoSquare(ensemblesToUse,:))
% imagesc(meanTSSquare(:,:))
title('stim')

colormap rdbu
caxis([-0.2 0.2])

ax(1)=subplot(2,2,4);
fillPlot(meanTSSquare(ensemblesToUse,:),[],'ci');

subplot(2,2,1)
imagesc(meanTSSquareNR(IndsUsed,:))
title('NoStim')

colormap rdbu
caxis([-0.1 0.1])
ax(2) = subplot(2,2,3);
fillPlot(meanTSSquareNR(IndsUsed,:),[],'ci');

linkaxes(ax);
