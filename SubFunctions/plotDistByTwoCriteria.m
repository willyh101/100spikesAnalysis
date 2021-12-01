function plotDistByTwoCriteria(All,outVars,opts,figNum)

plotTraces=opts.plotTraces;

%things to hold constant
ensemblesToUse = opts.ensemblesToPlot; %outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<0.5; % & outVars.meanEnsOSI<0.5 ;%  &  outVars.ensOSI>0.75;;

criteria = opts.criteria; %outVars.meanEnsOSI;% outVars.ensOSI;outVars.meanEnsOSI;
% useableCriteria = criteria(ensemblesToUse);
bins = opts.criteriaBins; %[0 0.5 inf];

criteria2 =   opts.criteria2; %outVars.ensOSI; %outVars.ensMaxD;
% useableCriteria = criteria2(ensemblesToUse);
% bins2 = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins2 = opts.criteria2Bins; %[0 0.3 0.7 inf];



numEns = numel(ensemblesToUse);

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
    
    if ensemblesToUse(i)
        ind = ensIndNumber(i);
        
        cellToUseVar = ~outVars.offTargetRiskEns{i} &  All(ind).out.anal.cellsToInclude;
        
        if opts.useVisCells
            cellToUseVar = cellToUseVar...
                & outVars.pVisR{ind} < 0.05;
        end
        if opts.useTunedCells
            cellToUseVar = cellToUseVar ...
                & outVars.osi{ind} > 0.25 ...
                ;
         end
        
        popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
        
    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')

figure(figNum);clf

colorListOri = colorMapPicker(numel(bins2),'plasma');
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
            [opts.criteriaName ': ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            [opts.criteria2Name ': ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            ['Num Ens: ' num2str(sum(ensembleSelecter )) ] ...
            } )
        
        ylabel('Pop Response (Mean \DeltaF/F)')
        axis square
        %     figure(11); tempax = subplot(1,1,1);
        %     [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
        %     eHandle{1}.Color = colorListOri{k};
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.25 0.3])

% p1 = ranksum(squeeze(dataForStats(3,1,:)),squeeze(dataForStats(3,3,:)));
% % [a p1] = ttest2(squeeze(dataForStats(3,1,:)),squeeze(dataForStats(3,3,:)));
% disp(['Co Tuned Close vs Far p = ' num2str(p1)]);
%
% p2 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(1,3,:)));
% % [a p2] = ttest2(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(1,3,:)));
% disp(['UnTuned Close vs Far p = ' num2str(p2)]);
%
% p3 = ranksum(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(3,1,:)));
% % [a p3] = ttest2(squeeze(dataForStats(1,1,:)),squeeze(dataForStats(3,1,:)));
% disp(['Close Tuned v Untuned p = ' num2str(p3)]);
%
% p4 = ranksum(squeeze(dataForStats(1,3,:)),squeeze(dataForStats(3,3,:)));
% % [a p4] = ttest2(squeeze(dataForStats(1,3,:)),squeeze(dataForStats(3,3,:)));
% disp(['Far Tuned v Untuned p = ' num2str(p4)]);
