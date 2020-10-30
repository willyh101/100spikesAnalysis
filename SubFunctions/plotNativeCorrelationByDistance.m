function plotNativeCorrelationByDistance(All,outVars)
%this function originally created to compare to Jenny Experiments. But it
%might be useful in the future. At present (Oct 8 2020), It doesn't involve
%any holography

indList = find(outVars.visPercent>0.1);


TCAll =[];
DiAll =[];
cellCounter=0;
for L = 1:numel(indList)
    ind = indList(L);
    disp(['Working on ind ' num2str(ind)]);

    cellsIncluded = All(ind).out.anal.pVisR<0.01;
%     cellToUse = cellsIncluded
smallResp = All(ind).out.anal.oriCurve(:,cellsIncluded);

pxToMu = 800/511;
    
    xyLoc = All(ind).out.exp.allCoM(cellsIncluded,:).*pxToMu;
    zLoc = (All(ind).out.exp.allDepth(cellsIncluded)-1)*30; %assumes 30um spacing
    
    xyzLoc = [xyLoc zLoc];
    
    Dist = nan([sum(cellsIncluded) sum(cellsIncluded)]);
    for i=1:sum(cellsIncluded);
        for k = 1:sum(cellsIncluded);
            if i~=k & i<k
                Dist(i,k) = sqrt( sum((xyzLoc(i,:) - xyzLoc(k,:)).^2));
            end
        end
    end
    
    tuneCorr = corr(smallResp);
    temp =eye(size(tuneCorr,1));
    tuneCorr(logical(temp(:)))=nan;
    
    TCAll = cat(1,TCAll,tuneCorr(:));
    DiAll = cat(1,DiAll,Dist(:));
    cellCounter = cellCounter  + sum(cellsIncluded);
    
end

figure(24);clf
    incr =25;
    distBins = [0:incr:900];
    
    corrByDistM = zeros([1 numel(distBins)-1]);
    corrByDistSE =  zeros([1 numel(distBins)-1]);
    for i=1:numel(distBins)-1
        corrToUse =TCAll(DiAll>distBins(i) & DiAll<=distBins(i+1));
        corrByDistM(i) = nanmean(corrToUse);
        corrByDistSE(i) = nanstd(corrToUse)/sqrt(numel(corrToUse));
        
    end
    e = errorbar(distBins(1:end-1)+incr/2,corrByDistM,corrByDistSE);
    e.Color = 'k';% CL{L};

ylim([-0 1])
xlim([0 750])
ylabel('Pairwise Signal Correlation (pearsons''s)')
xlabel('Pairwise Distance (\mum)')