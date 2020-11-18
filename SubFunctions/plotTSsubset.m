function [All outVars] = plotTSsubset(All,outVars)

%% Create time series plot

numExps = numel(All);
ensemblesToUse = outVars.ensemblesToUse; 
IndsUsed = outVars.IndsUsed;

minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));

clear allMeanTS allStdTS allnumTS allMeanTSVis allStdTSVis allnumTSVis
for ind=1:numExps
    us = unique(All(ind).out.exp.stimID);

    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    strtFrame = All(ind).out.anal.recStartFrame;
    newStart = strtFrame-minStrtFrame+1;
    
    pVisR = All(ind).out.anal.pVisR;
    osis = All(ind).out.anal.osi;
    
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
            cellsToUse = ~ROIinArtifact' & pVisR<0.05;
            dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end,trialsToUse &...
                All(ind).out.exp.stimID==s);
            mDat = mean(dat,3);
            mmDat = mean(mDat,1);
            meanVis = mmDat;
            sdDat = std(mDat);
            nDat = size(mDat,1);
            
        else
            exptEns = outVars.uExptsEns == ind;
            ensTuning = outVars.ensPO(exptEns);
            thisEnsTuning = ensTuning(i-1);
            
            cellsTuning = idx2ori(outVars.prefOris{ind}, [nan 0:45:315]);
%             coTunedCellsToEns = (mod(cellsTuning-90, 315) == thisEnsTuning);
            coTunedCellsToEns = (cellsTuning == thisEnsTuning);
        
            cellsToUse = ~ROIinArtifact' &...
                ~offTargetRisk(h,:) &...
                ~ismember(cellList,tg) &...
                pVisR<0.05 &...
                osis>0.5 &...
                coTunedCellsToEns;

            
            dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end,trialsToUse &...
                All(ind).out.exp.stimID==s);
            mDat = mean(dat,3);
            mmDat = mean(mDat,1)-meanVis;
            sdDat = std(mDat);
            nDat = size(mDat,1);
        end

        mRespTS(i,:) = mmDat; % mean response time series
        sRespTS(i,:) = sdDat; % std response time series (by cell);
        nResp(i) = nDat;
        
    end

    allMeanTS{ind} = mRespTS;
    allStdTS{ind} = sRespTS;
    allnumTS{ind} = nResp;
    
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

clim = [-0.4 0.4];

numCellsEachEns = outVars.numCellsEachEns;
uniqueNumCells = unique(numCellsEachEns(ensemblesToUse));
numSizes = numel(uniqueNumCells);

figure(2552);clf

clear ix ax
ix(1) =subplot(2,numSizes+1,1);
imagesc(meanTSSquareNR(IndsUsed,:))
title('Visual Response')
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
sgtitle('Pyramidal Cells Iso to Stimulus')

%co plot
figure(2992);clf


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