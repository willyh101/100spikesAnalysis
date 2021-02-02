function [All outVars] = createTSPlotAllVis(All,outVars)

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

    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    strtFrame = All(ind).out.anal.recStartFrame;
    newStart = strtFrame-minStrtFrame+1;
    
    pVisR = All(ind).out.anal.pVisR;
    
    trialsToUse = All(ind).out.anal.defaultTrialsToUse; 
    
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
            cellsToUse = ~ROIinArtifact' & pVisR<0.05 ;
        else
            cellsToUse = ~ROIinArtifact' &...
                ~offTargetRisk(h,:) &...
                ~ismember(cellList,tg) & pVisR<0.05;
        end
        
        dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end,trialsToUse &...
            All(ind).out.exp.stimID==s);

        mDat = mean(dat,3);
        mmDat = mean(mDat,1); %pop Average
        sdDat = std(mDat);
        nDat = size(mDat,1);


        mRespTS(i,:) = mmDat; % mean response time series
        sRespTS(i,:) = sdDat; % std response time series (by cell);
        nResp(i) = nDat;
        
    end
    All(ind).out.anal.mRespTSAllVis= mRespTS;
    All(ind).out.anal.sRespTSAllVis= sRespTS;
    All(ind).out.anal.nRespAllVis = nResp;
    
    
    allMeanTS{ind} = mRespTS;
    allStdTS{ind} = sRespTS;
    allnumTS{ind} = nResp;
    
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

clim = [-0.15 0.15];

figure(251);clf
ix(1) =subplot(2,2,2);
imagesc(meanTSSquare(ensemblesToUse,:))
ylabel('Ensemble')
xlabel('Frame')

% imagesc(meanTSSquare(:,:))
title('stim')

colormap rdbu
caxis(clim)

clear ax
ax(1)=subplot(2,2,4);

rc =rectangle('position',[minStrtFrame -0.025 6 0.045]);
rc.FaceColor = [rgb('FireBrick') 0.25];
rc.LineStyle = 'none';
hold on

fillPlot(meanTSSquare(ensemblesToUse,:),[],'ci');
ylabel('\DeltaZ-Scored dF/F')
xlabel('Frame')

% line([minStrtFrame minStrtFrame], [min(meanTSSquare(ensemblesToUse,:)) max(meanTSSquare(ensemblesToUse,:))]);
r = refline(0);
r.Color = rgb('grey');
r.LineStyle=':';
r.LineWidth=2;


ix(2) =subplot(2,2,1);
imagesc(meanTSSquareNR(IndsUsed,:))
title('NoStim')
ylabel('Mouse')
xlabel('Frame')

colormap rdbu
caxis(clim)
ax(2) = subplot(2,2,3);
if numel(IndsUsed)>1
fillPlot(meanTSSquareNR(IndsUsed,:),[],'ci');
else
    p= plot(meanTSSquareNR(IndsUsed,:));
    p.Color = 'k';
    p.LineWidth = 2;
end

r = refline(0);
r.Color = rgb('grey');
r.LineStyle=':';
r.LineWidth=2;
ylabel('\DeltaZ-Scored dF/F')
xlabel('Frame')
linkaxes(ax);
linkaxes([ix ax],'x')
xlim([1 size(meanTSSquareNR,2)]);