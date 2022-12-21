function [popVal] = subsetPopResponseMatchedNoStim(All,outVars,opts)

ensemblesToUse = opts.ensemblesToPlot; 

numEns = numel(outVars.ensIndNumber);

if ~isfield(opts,'visCond') || opts.visCond==1;
    skipVis=1;
    v=1;
else
    skipVis = 0;
    v = opts.visCond;
end

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end
    
    if ensemblesToUse(i)
            ind = outVars.ensIndNumber(i);
        
        cellToUseVar = ~outVars.offTargetRiskEns{i} &  All(ind).out.anal.cellsToInclude & ~All(ind).out.anal.ROIinArtifact';
        
        if opts.useVisCells
            cellToUseVar = cellToUseVar...
                & outVars.pVisR{ind} < 0.05;
        end
        if opts.useTunedCells
            cellToUseVar = cellToUseVar ...
                & outVars.osi{ind} > 0.25 ...
                ;
        end
        
        if isfield(opts,'variableCellFun') & ~isempty(opts.variableCellFun)
            cellToUseVar = cellToUseVar ...
                & eval(opts.variableCellFun) ;
        end
        
        if sum(cellToUseVar)<opts.minNumberOfCellsPerCondition
            cellToUseVar = logical(zeros(size(cellToUseVar)));
        end

         
            s = 1; %outVars.ensHNumber(i);
            
            respMat = All(ind).out.anal.respMat;
            baseMat = All(ind).out.anal.baseMat;
            
            if skipVis
                popVal(i) = nanmean(squeeze(respMat(s,1,cellToUseVar) - baseMat(s,1,cellToUseVar)));
            else
                popVal(i) = nanmean(squeeze(respMat(s,v,cellToUseVar) - baseMat(s,v,cellToUseVar)));
            end
            
%         popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popVal(i) = nan;
    end
end
disp('Done')