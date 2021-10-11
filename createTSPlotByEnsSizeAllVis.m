function [All outVars] = createTSPlotByEnsSizeAllVis(All,outVars)

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
    
    if isfield(All(ind).out.exp, 'dataToUse')
        dataToUse = All(ind).out.exp.dataToUse;
    else
        disp(['ind ' num2str(ind) '. no data to use, using zdfData']);
        dataToUse = All(ind).out.exp.zdfData;
    end
    
    
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
        
        if  h == 0
            cellsToUse = ~ROIinArtifact' & pVisR<0.05 ;
        else
            cellsToUse = ~ROIinArtifact' &...
                ~offTargetRisk(h,:) &...
                ~ismember(cellList,tg) &...
                pVisR<0.05;
        end

        dat = dataToUse(cellsToUse,newStart:end,trialsToUse &...
            All(ind).out.exp.stimID==s);

        mDat = mean(dat,3);
        mmDat = mean(mDat,1); %pop Average
        sdDat = std(mDat);
        nDat = size(mDat,1);


        mRespTS(i,:) = mmDat; % mean response time series
        sRespTS(i,:) = sdDat; % std response time series (by cell);
        nResp(i) = nDat;
        
    end


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

clim = [-0.3 0.3];

numCellsEachEns = outVars.numCellsEachEns;
uniqueNumCells = unique(numCellsEachEns(ensemblesToUse));
numSizes = numel(uniqueNumCells);

figure(252);clf

clear ix ax
ix(1) =subplot(2,numSizes+1,1);
imagesc(meanTSSquareNR(IndsUsed,:))
title('NoStim')
ylabel('Mouse')
xlabel('Frame')

colormap rdbu
caxis(clim)

ax(1) = subplot(2,numSizes+1,numSizes+2);
if numel(IndsUsed)>1
    lineCol = rgb('dimgrey');
    faceCol = rgb('dimgrey');
    fAlpha = 0.5;
    
    fillPlot(meanTSSquareNR(IndsUsed,:),[],'ci',...
        lineCol,...
        'none',...
        faceCol,...
        fAlpha);
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

colorList = colorMapPicker(numSizes,outVars.defaultColorMap);

for i = 1:numSizes
    ix(end+1) =subplot(2,numSizes+1,1+i);
    imagesc(meanTSSquare(ensemblesToUse & numCellsEachEns==uniqueNumCells(i),:))
    ylabel('Ensemble')
    xlabel('Frame')
    
    % imagesc(meanTSSquare(:,:))
    title([num2str(uniqueNumCells(i)) ' Cell Stim'])
    
    colormap rdbu
    caxis(clim)

    ax(end+1)=subplot(2,numSizes+1,numSizes+2+i);
    
    rc =rectangle('position',[minStrtFrame -0.025 6 0.045]);
    rc.FaceColor = [rgb('FireBrick') 0.25];
    rc.LineStyle = 'none';
    hold on
    
    lineCol = colorList{i};
    faceCol = colorList{i};
    fAlpha = 0.5;
    
    fillPlot(meanTSSquare(ensemblesToUse & numCellsEachEns==uniqueNumCells(i),:),[],'ci',...
        lineCol,...
        'none',...
        faceCol,...
        fAlpha);
    ylabel('\DeltaZ-Scored dF/F')
    xlabel('Frame')
    
    % line([minStrtFrame minStrtFrame], [min(meanTSSquare(ensemblesToUse,:)) max(meanTSSquare(ensemblesToUse,:))]);
    r = refline(0);
    r.Color = rgb('grey');
    r.LineStyle=':';
    r.LineWidth=2;
end
linkaxes([ix ax],'x')
linkaxes(ax);
xlim([1 size(meanTSSquareNR,2)]);

%co plot
figure(292);clf


if numel(IndsUsed)>1
    lineCol = rgb('dimgrey');
    faceCol = rgb('dimgrey');
    fAlpha = 0.5;
    
    fillPlot(meanTSSquareNR(IndsUsed,:),[],'ci',...
        lineCol,...
        'none',...
        faceCol,...
        fAlpha);
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

rc =rectangle('position',[minStrtFrame -0.025 6 0.045]);
    rc.FaceColor = [rgb('FireBrick') 0.25];
    rc.LineStyle = 'none';
    hold on

for i = 1:numSizes


    
    lineCol = colorList{i};
    faceCol = colorList{i};
    fAlpha = 0.5;
    
    fillPlot(meanTSSquare(ensemblesToUse & numCellsEachEns==uniqueNumCells(i),:),[],'ci',...
        lineCol,...
        'none',...
        faceCol,...
        fAlpha);
    ylabel('\DeltaZ-Scored dF/F')
    xlabel('Frame')
    
    % line([minStrtFrame minStrtFrame], [min(meanTSSquare(ensemblesToUse,:)) max(meanTSSquare(ensemblesToUse,:))]);
    r = refline(0);
    r.Color = rgb('grey');
    r.LineStyle=':';
    r.LineWidth=2;
end