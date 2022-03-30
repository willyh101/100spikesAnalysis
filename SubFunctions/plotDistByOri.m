function plotDistByOri(All,outVars,opts)
% opts.distBins =[0:25:150];% [0:25:400];
% 
% opts.distType = 'harm'; 
% opts.distBins = [0:50:400];
names = outVars.names;

[outVars] = grandDistanceMaker(opts,All,outVars);
numEns = numel(outVars.ensStimScore);

plotOrientation = opts.plotOrientation; %as opposed to Direction

minNumberOfCellsPerCondition = opts.minNumberOfCellsPerCondition; %set to -1 to ignore

%this is where you change the criteria of what ensembles are included
ensemblesToUse = opts.ensemblesToPlot; %outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI>0.7;% & outVars.meanEnsOSI>0.25;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse &  outVars.ensOSI>0.7 & outVars.meanEnsOSI>0.5;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;
% ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 &  outVars.ensOSI>0.7;% & outVars.meanEnsOSI>0.25;% & outVars.numMatchedTargets>=3 & outVars.ensMaxD>-475;% outVars.numMatchedTargets>=3 &    ;


oriVals = [NaN 0:45:315];
% numEns = numel(outVars.posCellbyInd);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

% disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
% disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
% disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);

clear popToPlot cellsPerEnsCount
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
% 
    diffsPossible = [0 45 90 135 180];
    
    if ensemblesToUse(i)
        ind = ensIndNumber(i);
        
        cellOris = oriVals(outVars.prefOris{ind});
        cellOrisDiff = abs(cellOris-outVars.ensPO(i));
        cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
        %
        if plotOrientation
            cellOrisDiff(cellOrisDiff==135)=45;
            cellOrisDiff(cellOrisDiff==180)=0;
            diffsPossible = [0 45 90];
        end
        
        %%This is where you change the criteria of what cells are included
        cellToUseVar = ~outVars.offTargetRiskEns{i}...
            & outVars.pVisR{ind} < 0.05 ...
            & outVars.osi{ind} > 0.25 ...
            & All(ind).out.anal.cellsToInclude ...
            ...& outVars.posCellbyInd{i} ... %if you want to only include cells that went up
            ...& outVars.isRedByEns{i} ...  %if you want to excluded red cells (i.e. interneurons)
        ;

    for k=1:numel(diffsPossible)
        
        cells2Plot =cellToUseVar & abs(cellOrisDiff)==diffsPossible(k);
        cellsPerEnsCount(i,k) = sum(cells2Plot);
    end
    for k=1:numel(diffsPossible)
        
        if all(cellsPerEnsCount(i,:)>minNumberOfCellsPerCondition)
            cells2Plot =cellToUseVar & abs(cellOrisDiff)==diffsPossible(k);
            popToPlot(i,:,k) = popDistMakerSingle(opts,All(ensIndNumber(i)),cells2Plot,0,ensHNumber(i));
        else 
            cellsPerEnsCount(i,:) =zeros([1 numel(diffsPossible)]);
        end        
    end
    %plot all nonVis cells
    cellToUseVar = ~outVars.offTargetRiskEns{i}...
        & outVars.pVisR{ind} > 0.05 ...
        & All(ind).out.anal.cellsToInclude ...
        ;
    cells2Plot =cellToUseVar;
    popToPlot(i,:,k+1) = popDistMakerSingle(opts,All(ensIndNumber(i)),cells2Plot,0,ensHNumber(i));
            
    else
        popToPlot(i,:,:) = nan([numel(opts.distBins)-1 1 numel(diffsPossible)+1]);
    end
end
disp('Done')
sum(cellsPerEnsCount(:,1)>0);
ensUsedActual = cellsPerEnsCount(:,1)>0;


disp(['Total of ' num2str(sum(ensUsedActual)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensUsedActual)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensUsedActual)))))) ' Mice']);


%%Plot Dist REsp
figure(10);clf
opts.distAxisRange = [0 350];

figure(11);clf;hold on
figure(14);clf

% ax =subplot(1,1,1);
%  hold on
colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
colorListOri{end+1} = rgb('grey');
dataForStats=[];
clear ax

if plotOrientation
    diffsPossible = [0 45 90];
end

plotNonVis=0;
if plotNonVis
numToPlot = numel(diffsPossible)+1;
else
    numToPlot = numel(diffsPossible);
end

for k = 1:numToPlot
    figure(10);
    ax(k) =subplot(1,numToPlot,k);
    if k<=numel(diffsPossible)
    title(['Cells Pref Angle \Delta' num2str(diffsPossible(k)) '\circ'])
    else
            title(['Non Visually Responsive Cells'])

    end
    [eHandle outData] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,ax(k));
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end
    eHandle{1}.CapSize =0;
    ylabel('Pop Response (Mean \DeltaF/F)')

            dataForStats(k,:) = outData{1}.dat(:,1);

    
    figure(11); tempax = subplot(1,1,1);
    [eHandle] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,tempax);
    eHandle{1}.Color = colorListOri{k};
    ylabel('Pop Response (Mean \DeltaF/F)')

    if k ==1 || k==3
        figure(14); tempax = subplot(1,1,1);
        [eHandle] = plotDistRespGeneric(popToPlot(:,:,k),outVars,opts,tempax);
        eHandle{1}.Color = colorListOri{k};
        eHandle{1}.CapSize =0;
        ylabel('Pop Response (Mean \DeltaF/F)')

        if k==1
        delete(eHandle{end})
        end
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])

figure(10);
xlim([0 opts.distBins(end)])
ylim([-0.25 0.25])

figure(14)
xlim([0 opts.distBins(end)])
legend('Iso','Ortho')



figure(12);clf
datToPlot = squeeze(popToPlot(ensemblesToUse,1,1:numToPlot));
ciToPlot = nanstd(datToPlot,[],1)./sqrt(size(datToPlot,1))*1.93;
errorbar(nanmean(datToPlot),ciToPlot);
p = anova1(datToPlot,[],'off');
% p2 = signrank(datToPlot(:,1),datToPlot(:,3));
% p2 = signrank([datToPlot(:,1);datToPlot(:,5)],datToPlot(:,3));
% [~,p2] = ttest2(datToPlot(:,1),datToPlot(:,3));
% [~,p2] = ttest2([datToPlot(:,1);datToPlot(:,5)],datToPlot(:,3));

[~, p2T] = ttest2(squeeze(dataForStats(1,:)),squeeze(dataForStats(3,:)));
disp(['ranksum iso v ortho Ttest p = ' num2str(p2T)]);


p2 = ranksum(squeeze(dataForStats(1,:)),squeeze(dataForStats(3,:)));
disp(['ranksum iso v ortho p = ' num2str(p2)]);


title({['Pvalue ' num2str(p)]...
    ['Iso vs Ortho PValue: ' num2str(p2) ' Ttest: ' num2str(p2T)] })

% title('Cells by Tuning')
% ax2 =subplot(1,2,2);
% plotDistRespGeneric(popToPlotNeg,outVars,opts,ax2);
% title('Cells That go down')
