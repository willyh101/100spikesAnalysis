function plotSpaceAndFeature(All,outVars,opts,figNum)
numEns = numel(outVars.ensStimScore);

ensemblesToUse = ...
    outVars.ensemblesToUse ...
    & outVars.numCellsEachEnsBackup==10 ...
    ...& outVars.meanEnsOSI>0.5 ...
    ;
criteria =  outVars.ensMaxD; outVars.ensMeaD;%   outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse);
% bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins = [0 400 500 inf];


% bins = [linspace(min(useableCriteria),max(useableCriteria),5)];
criteria2 =   outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria2(ensemblesToUse);
% bins2 = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins2 = [0 0.3 0.7 inf];

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
    
    if ensemblesToUse(i)
        ind = ensIndNumber(i);
        
        cellToUseVar = ~outVars.offTargetRiskEns{i}...
              & outVars.pVisR{ind} < 0.05 ...
              & All(ind).out.anal.cellsToInclude ...
            ... & All(ind).out.red.isRed ...
            ... & outVars.osi{ind} > 0.25 ...
            ... & outVars.posCellbyInd{i} ...
            ;
        sum(cellToUseVar);
        
        popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
        
    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')

figure(figNum);clf

colorListOri = colorMapPicker(numel(bins),'plasma');
clear ax
dataForStats =[];
c=0;
for i = 1:numel(bins2)-1
    for k = 1:numel(bins)-1
        c=c+1;
        
        ax(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
        popToPlotTemp = popToPlot;
        ensembleSelecter  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;
        popToPlotTemp(~ensembleSelecter,:)=NaN;
        [eHandle outDat] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(c));
        eHandle{1}.CapSize =0;
        if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
            eHandle{1}.Color = colorListOri{k};
        end
        dataForStats(i,k,:) = outDat{1}.dat(:,1);
        
        if opts.plotTraces
            hold on
            distBins = opts.distBins;
            distBinSize = distBins(2)-distBins(1);
            plot(distBins(2:end)-distBinSize/2,outDat{1}.dat','color',rgb('grey'));
        end
        
        title({...
            ['Ens Dist: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            ['Ens OSI: ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            ['Num Ens: ' num2str(sum(ensembleSelecter )) ] ...
            } )
        
        ylabel('Pop Response (Mean \DeltaF/F)')
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.25 0.3])
try
p1 = ranksum(squeeze(dataForStats(3,1,:)),squeeze(dataForStats(3,3,:)));
disp(['Co Tuned Close vs Far p = ' num2str(p1)]);
p2 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(1,3,:)));
disp(['UnTuned Close vs Far p = ' num2str(p2)]);

p3 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(3,1,:)));
disp(['Close Tuned v Untuned p = ' num2str(p3)]);
p4 = ranksum(squeeze(dataForStats(1,3,:)),squeeze(dataForStats(3,3,:)));
disp(['Far Tuned v Untuned p = ' num2str(p4)]);
catch
end