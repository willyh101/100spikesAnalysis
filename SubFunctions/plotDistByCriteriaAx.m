function [eHandle outDat] = plotDistByCriteriaAx(All,outVars,opts,ax)

plotTraces=opts.plotTraces;

%things to hold constant
ensemblesToUse = opts.ensemblesToPlot; %outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<0.5; % & outVars.meanEnsOSI<0.5 ;%  &  outVars.ensOSI>0.75;;
criteria = opts.criteriaToSplit; %outVars.meanEnsOSI;% outVars.ensOSI;outVars.meanEnsOSI;
% bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
% bins = [0 prctile(useableCriteria,33) prctile(useableCriteria,66) max(useableCriteria)];

% bins = [0 0.25 0.5 0.75 max(useableCriteria)];
% bins = [0 0.33 0.7 max(useableCriteria)];
% bins = [0 0.3 0.7 max(useableCriteria)];
bins = opts.criteriaBins; %[0 0.5 inf];

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
        
        cellToUseVar = ~outVars.offTargetRiskEns{i} &  All(ind).out.anal.cellsToInclude & ~All(ind).out.anal.ROIinArtifact';
        
        if ~isempty(opts.variableCellFun)
            try
            cellToUseVar = cellToUseVar & eval(opts.variableCellFun);
            catch
                disp(['Variable Cell Fun Error. ens: ' num2str(i)])
                cellToUseVar= zeros(size(cellToUseVar));
            end
        end
        
        if opts.useVisCells == 1
            cellToUseVar = cellToUseVar...
                & outVars.pVisR{ind} < 0.05;
        elseif opts.useVisCells == -1
            cellToUseVar = cellToUseVar...
                & outVars.pVisR{ind} > 0.1;
        end
        if opts.useTunedCells
            cellToUseVar = cellToUseVar ...
                & outVars.osi{ind} > 0.25 ...
                ;
        end
        
        
        if all(cellToUseVar==0)
            popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
        else
            popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
        end
    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')


colorListOri = colorMapPicker(numel(bins)-1,'plasma');
for k = 1:numel(bins)-1
  
%     figure(figNum+1);
    ensembleSelecter = criteria<bins(k) | criteria>bins(k+1)  | ~ensemblesToUse | isnan(criteria);
    title({...
        ['Ens Criteria: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
        ['Num Ens: ' num2str(sum(~ensembleSelecter )) ] ...
        } );
    
    popToPlotTemp = popToPlot;
    popToPlotTemp(ensembleSelecter,:)=NaN;
    [eHandle outDat] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax);
    eHandle{1}.CapSize=0;
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end
    if plotTraces
        hold on
        distBins = opts.distBins;
        distBinSize = distBins(2)-distBins(1);
        plot(distBins(2:end)-distBinSize/2,outDat{1}.dat','color',rgb('grey'));
    end
    
end
% linkaxes(ax)
% figure(figNum+1);
% xlim([0 opts.distBins(end)])
% ylim([-0.1 0.3])
% figure(figNum+1);
% xlim([0 opts.distBins(end)])
% ylim([-0.1 0.3])