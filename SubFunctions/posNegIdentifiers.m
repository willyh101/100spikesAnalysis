function [All outVars] = posNegIdentifiers(All,outVars,opts)
if ~isfield(opts,'posNegThreshold')
    Threshold =0;
else
    Threshold = opts.posNegThreshold ;
end

numExps = numel(All);

posCellbyInd=[];
negCellbyInd=[];
for ind=1:numExps
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    numEns = size(respMat,1);
    clear posCells negCells
    for i=1:numEns;
    posCells(i,:)= (respMat(i,1,:)-baseMat(i,1,:)) >Threshold;
    negCells(i,:)= (respMat(i,1,:)-baseMat(i,1,:)) <Threshold*-1;
    
    if i>1
       posCellbyInd{end+1} = squeeze( posCells(i,:));
       negCellbyInd{end+1} = squeeze(negCells(i,:));
    end
    end
    
    All(ind).out.anal.posCells=posCells;
    All(ind).out.anal.negCells=negCells;
end
    
outVars.posCellbyInd=posCellbyInd;
outVars.negCellbyInd=negCellbyInd;
    

    