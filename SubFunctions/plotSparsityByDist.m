function [outData] = plotSparsityByDist(All,outVars,opts,ax)

subtractBaseline = opts.subtractBaseline;

ensemblesToUseList = find(opts.ensemblesToPlot );

for i=1:numel(ensemblesToUseList)
    if mod(i,50)==0
        disp('.')
    else
        fprintf('.')
    end
    
    ens = ensemblesToUseList(i);
    ind = outVars.ensIndNumber(ens);
    hNum = outVars.ensHNumber(ens);
    
    holo = All(ind).out.exp.stimParams.roi{hNum}; % Better Identifying ensemble
    Tg=All(ind).out.exp.rois{holo};
    
    
    if opts.useVisCells
        cellToUse = ~outVars.offTargetRiskEns{ens}...
            & outVars.pVisR{ind} < 0.05 ...
            ...& outVars.osi{ind} > 0.25 ...
            & All(ind).out.anal.cellsToInclude ...
            ... & outVars.isRedByEns{i} ...
            ;
    else
        cellToUse = ~outVars.offTargetRiskEns{ens}...
            & outVars.pVisR{ind} > 0.05 ...
            ...& outVars.osi{ind} > 0.25 ...
            & All(ind).out.anal.cellsToInclude ...
            ... & outVars.isRedByEns{i} ...
            ;
    end
   
    
    %the not by distance part
    theseData = squeeze(All(ind).out.anal.respMat(hNum,1,cellToUse));
    if subtractBaseline
        baseData = squeeze(All(ind).out.anal.baseMat(hNum,1,cellToUse));
        theseData = theseData-baseData;
    end
   
    L2 = sqrt(sum(theseData.^2))/numel(theseData);
    L1 = sum(abs(theseData))/numel(theseData);
    
    allL1(i) = L1;
    allL2(i) = L2;
    allSparse(i) = L2/L1;
    
    
    %by distance
    %First create stimDistance
    if ~isfield(All(ind).out.anal,'StimDistance')
        StimDistance = zeros([size(stimCoM,1) numCells]);
        for ii=1:size(stimCoM,1);
            for k=1:numCells;
                StimDistance(ii,k) = sqrt(sum((stimLoc(ii,:)-allLoc(k,:)).^2));
            end
        end
    else
        StimDistance = All(ind).out.anal.StimDistance;
    end
    
    dists = StimDistance(Tg,:);
    
    switch opts.distType
        case 'min'
            minDist = min(dists,[],1);
            distToUse = minDist;
        case 'geo'
            try
                geoDist = geo_mean(dists,1);
            catch
                geoDist = geomean(dists,1);
            end
            distToUse = geoDist;
        case 'mean'
            meanDist = mean(dists,1);
            distToUse = meanDist;
        case 'harm'
            harmDist = harmmean(dists,1);
            distToUse = harmDist;
        case 'median'
            medDist = median(dists,1);
            distToUse = medDist;
        case 'centroid'
            stimCentroid = mean(stimLoc(Tg,:));
            tempDists =[];
            for idx =1:numCells
                tempDists(idx) = sqrt(sum((stimCentroid-allLoc(idx,:)).^2));
            end
            distToUse = tempDists;
        otherwise
            disp('dist type not understood, using min')
            minDist = min(dists,[],1);
            distToUse =minDist;
    end
    
   roisInOrder =  [All(ind).out.exp.stimParams.roi{:}]';
   noStimIdx = find(roisInOrder==0);
    
    clear L2byDist L1byDist allSparse mnByDist
    for k = 1:numel(opts.distBins)-1
        d0 = opts.distBins(k);
        d1 = opts.distBins(k+1);
        
        cell2use = cellToUse & distToUse>=d0 & distToUse<d1;
        theseData = squeeze(All(ind).out.anal.respMat(hNum,1,cell2use));
        if subtractBaseline
            baseData = squeeze(All(ind).out.anal.baseMat(hNum,1,cell2use));
            theseData = theseData-baseData;
        end
        
        theseNoStimData = squeeze(All(ind).out.anal.respMat(noStimIdx,1,cell2use));
        if subtractBaseline
            baseData = squeeze(All(ind).out.anal.baseMat(noStimIdx,1,cell2use));
            theseNoStimData = theseNoStimData-baseData;

        end
%                                 theseNoStimData = baseData;

        %         theseData = theseData>0;
        
        L2 = sqrt(sum(theseData.^2))/numel(theseData);
        L1 = sum(abs(theseData))/numel(theseData);
        mnDat = nanmean(theseData);

        %
        %
        % %Rolls and Tovee
        % L1 = (sum(abs(theseData))/numel(theseData))^2;
        % L2 = sum(theseData.^2)/numel(theseData);
        
        L2NS = sqrt(sum(theseNoStimData.^2))/numel(theseNoStimData);
        L1NS = sum(abs(theseNoStimData))/numel(theseNoStimData);
        mnDatNS = nanmean(theseNoStimData);

        
        if sum(cell2use)==0
            L1byDist(k) = nan;
            L2byDist(k) = nan;
            allSparse(k) = nan;
            mnByDist(k) = nan;
            allSparseNS(k) = nan;

        else
            
            L1byDist(k) = L1;
            L2byDist(k) = L2;
            allSparse(k) = L2/L1;
            %             allSparse(k) = L1/L2;
            allSparseNS(k) = L2NS/L1NS;

            
            mnByDist(k) = mnDat;
            mnNSByDist(k) = mnDatNS;
        end
    end
    
    allL1byDist(:,i) = L1byDist;
    allL2byDist(:,i) = L1byDist;
    allSparseByDist(:,i) = allSparse;
    allMnByDist(:,i) = mnByDist;
    allSparseNSByDist(:,i) = allSparseNS;
        allMnNSByDist(:,i) = mnNSByDist;

end

outData.allL1 = allL1;
outData.allL2 = allL2;
outData.allSparse = allSparse;
outData.alllL1byDist = allL1byDist;
outData.alllL2byDist = allL2byDist;
outData.alllSparseByDist = allSparseByDist;
outData.alllSparseNSByDist = allSparseNSByDist;


if nargin>3
    subplot(ax);
    toPlot =allSparseByDist;% allMnByDist; ;%allSparseByDist;
    mn = nanmean(toPlot,2);
    sd = nanstd(toPlot,[],2);
    nm = sum(~isnan(toPlot),2);
    se = sd./sqrt(nm);
    
    e = errorbar(opts.distBins(1:numel(mn)),mn,se);
    e.CapSize=0;
    e.LineWidth=2;
    outData.figHandle = e;
    
    hold on
    toPlot =allSparseNSByDist;  allMnNSByDist; ;
    mn = nanmean(toPlot,2);
    sd = nanstd(toPlot,[],2);
    nm = sum(~isnan(toPlot),2);
    se = sd./sqrt(nm);
    
    e2 = errorbar(opts.distBins(1:numel(mn)),mn,se);
    e2.CapSize=0;
    e2.LineWidth=2;
    e2.Color = rgb('grey');
    
end



