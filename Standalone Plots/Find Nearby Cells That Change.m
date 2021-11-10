
 c=1;cellPairOfInterest=[]; cellPairOfInterestEnsemble=[];
for ind = IndsUsed
    
    eligibleEns = find(outVars.ensemblesToUse & ensIndNumber==ind);
    
    HTs=[];
    for i = 1:numel(eligibleEns)
        e = eligibleEns(i);
        s = outVars.ensHNumber(e); %i think this is the order of resp
        h = All(ind).out.exp.stimParams.roi{s};
        htg = All(ind).out.exp.holoTargets{h};
        %        htg(isnan(htg))=[];
        HTs(i,:) = htg;
    end
    
    uHTs = unique(HTs);
    uHTs(isnan(uHTs))=[];
    
    countHTs = [];
    for i=1:numel(uHTs)
        countHTs(i) = sum(HTs(:)==uHTs(i));
    end
    
    targsHitMulti = uHTs(countHTs>1);
    
    allDist = All(ind).out.anal.allDistance;
    range =25; %distance away to consider
    adjacentCells =[];
    for i=1:numel(targsHitMulti)
        t = targsHitMulti(i);
        adjacentCells = find(allDist(t,:)<=range & allDist(t,:)>0);
        adjacentCellsAll{i} = adjacentCells;
        
        es = eligibleEns(any(HTs==t,2));
        
        if ~isempty(adjacentCells)
            for k=1:numel(adjacentCells)
                a = adjacentCells(k);
                AdjCellResp=[];
                for j = 1:numel(es)
                    e = es(j);
                    s = outVars.ensHNumber(e);
                    h = All(ind).out.exp.stimParams.roi{s};
                    htg = All(ind).out.exp.holoTargets{h};
                    if ismember(a,htg)
                        AdjCellResp(j) = NaN;
                    else
                        AdjCellResp(j)=squeeze(All(ind).out.anal.respMat(s,1,a)-All(ind).out.anal.baseMat(s,1,a));
                    end
                end
                %                                AdjCellResp
                upperBound = 0.2;
                if min(AdjCellResp)<0 && max(AdjCellResp)>upperBound;
                    cellPairOfInterest(c,:) = [ind t a];
                    cellPairOfInterestEnsemble{c} = [es(AdjCellResp<0) es(AdjCellResp>upperBound)];
                    c=c+1;
                end
            end
        end
    end
    
    
end

cellPairOfInterest

%%
figure(10001);clf
 for i =1:size(cellPairOfInterest,1)
     clf
     ind = cellPairOfInterest(i,1);
     t = cellPairOfInterest(i,2);
     a = cellPairOfInterest(i,3);
     
     e1 = cellPairOfInterestEnsemble{i}(1);
     e2 = cellPairOfInterestEnsemble{i}(2);

     try
     ax = subplot(2,6,[1 7]);
     [All(ind)] = plotCellMask(All(ind),t,ax);
     title(num2str(t))
     catch;end
     
     ax = subplot(2,6,2);
     plotCaRaster(All,outVars,e1,t,ax)
     ax = subplot(2,6,8);
     plotCaRaster(All,outVars,e2,t,ax)
     
     ax = subplot(2,6,3);
     plotCellMean(All,outVars,e1,t,rgb('darkmagenta'),ax)
     ax = subplot(2,6,9);
     plotCellMean(All,outVars,e2,t,rgb('darkmagenta'),ax)

     try
     ax = subplot(2,6,[4 10]);
     [All(ind)] = plotCellMask(All(ind),a,ax);
         title(num2str(a))
     catch; end
     
         ax = subplot(2,6,5);
     plotCaRaster(All,outVars,e1,a,ax)
     ax = subplot(2,6,11);
     plotCaRaster(All,outVars,e2,a,ax)

          ax = subplot(2,6,6);
     plotCellMean(All,outVars,e1,a,rgb('darkmagenta'),ax)
     ax = subplot(2,6,12);
     plotCellMean(All,outVars,e2,a,rgb('darkmagenta'),ax)
     
     pause
 end
 
 disp('done')