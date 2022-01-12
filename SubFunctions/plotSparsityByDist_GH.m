function [outData] = plotSparsityByDist_GH(All,outVars,opts,ax)

subtractBaseline = opts.subtractBaseline;

ensemblesToUseList = find(opts.ensemblesToPlot );

omitFlag =0;


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
    
%     try
%         Tg=All(ind).out.exp.rois{holo+1};
%     catch
%         Tg=All(ind).out.exp.rois{holo-1};
%     end

% try rand targets
% Tg = randi(numel(All(ind).out.exp.targetedCells),[1 10]);
    
% Tg = randi(numel(outVars.pVisR{ind}),[1 10]);
     
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
    
    numCells = All(ind).out.anal.numCells;
    allLoc = [All(ind).out.exp.allCoM (All(ind).out.exp.allDepth-1)*30];
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
    
%     allCoM = All(ind).out.exp.allCoM;
%     StimDistance = zeros([size(allCoM,1) numCells]);
%     for ii=1:size(allCoM,1);
%         for k=1:numCells;
%             StimDistance(ii,k) = sqrt(sum((allLoc(ii,:)-allLoc(k,:)).^2));
%         end
%     end
    
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
        if opts.subSampleN>0
            if sum(cell2use) < opts.subSampleN
                %                 disp('Error Not enough Cells this bin.. blanking')
                omitFlag = omitFlag+1;
                cell2use = logical(zeros(size(cell2use)));
            else
                cellList = find(cell2use);
                newCells = cellList(randperm(numel(cellList),opts.subSampleN));
                cell2use = zeros(size(cell2use));
                cell2use(newCells)=1;
                cell2use = logical(cell2use);
                
            end
        end
        
        if opts.minSampleN>0
            if sum(cell2use) < opts.minSampleN
                omitFlag = omitFlag+1;
                cell2use = logical(zeros(size(cell2use)));
            end
        end
        %         disp(num2str(sum(cell2use)));
        
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
        
        
        switch opts.sparseAlgo
            case 'L2/L1'
                L2 = sqrt(sum(theseData.^2))/numel(theseData);
                L1 = sum(abs(theseData))/numel(theseData);
                
                L2NS = sqrt(sum(theseNoStimData.^2))/numel(theseNoStimData);
                L1NS = sum(abs(theseNoStimData))/numel(theseNoStimData);
                
                if sum(cell2use)==0
                    allSparse(k) = nan;
                    allSparseNS(k) = nan;
                else
                    allSparse(k) = L2/L1;
                    allSparseNS(k) = L2NS/L1NS;
                end
            case 'popKurtosis'
                
                num_to_sample = 20;
                if numel(theseData)>=num_to_sample
                    theseData = theseData(1:num_to_sample);
                    
                    mn = mean(theseData);
                    sd = std(theseData);
                    nm = numel(theseData);
                    
                    Kp = sum( ((theseData-mn)/sd) .^4)/nm - 3;
                    %                 Kp = (sum((theseData-mn)/sd ).^4 )/nm - 3;
                    
                    allSparse(k) = Kp;                  
                else
                    allSparse(k) = nan;
                end
                
                if numel(theseNoStimData)>=num_to_sample
                    
                   theseNoStimData = theseNoStimData(1:num_to_sample);
                    
                    mn = mean(theseNoStimData);
                    sd = std(theseNoStimData);
                    nm = numel(theseNoStimData);
                    
                    KpNS = sum( ((theseNoStimData-mn)/sd) .^4)/nm - 3;
                    %                 KpNS = (sum((theseNoStimData-mn)/sd ).^4 )/nm - 3;
                    
                    allSparseNS(k) =KpNS;
                else               
                    allSparseNS(k) = nan;
                end
                
                if sum(cell2use)==0
                    allSparse(k) = nan;
                    allSparseNS(k) = nan;
                end

            case 'treves-Rolls'
                
                if numel(theseData)>10
                    theseData = theseData(1:10);
                    
                    mn = mean(theseData);
                    sd = std(theseData);
                    nm = numel(theseData);
                    
                    S = 1 - (( sum(abs(theseData)./nm)^2 ) / (sum(theseData.^2 / nm)) );
                    
                    mn = mean(theseNoStimData);
                    sd = std(theseNoStimData);
                    nm = numel(theseNoStimData);
                    
                    SNS = 1 - (( sum(abs(theseNoStimData)./nm)^2 ) / (sum(theseNoStimData.^2 / nm)) );
                    
                    if sum(cell2use)==0
                        allSparse(k) = nan;
                        allSparseNS(k) = nan;
                    else
                        allSparse(k) = S;
                        allSparseNS(k) =SNS;
                    end
                else
                    allSparse(k) = nan;
                    allSparseNS(k) = nan;
                end
                
                
            case 'mean'
                if sum(cell2use)==0
                    allSparse(k) = nan;
                    allSparseNS(k) = nan;
                else
                    allSparse(k) = mean(theseData);
                    allSparseNS(k) =mean(theseNoStimData);
                end
            otherwise
                disp('did not understand sparse algo')
        end
    end
    
    
    
    %the inportant ones
    allSparseByDist(:,i) = allSparse;
    allSparseNSByDist(:,i) = allSparseNS; %no stim
    
end

if opts.subSampleN>0 || opts.minSampleN>0
    disp([num2str(omitFlag) ' bin(s) did not have enough data to include'])
end



outData.alllSparseByDist = allSparseByDist;


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
    toPlot =allSparseNSByDist;
    mn = nanmean(toPlot,2);
    sd = nanstd(toPlot,[],2);
    nm = sum(~isnan(toPlot),2);
    se = sd./sqrt(nm);
    
    e2 = errorbar(opts.distBins(1:numel(mn)),mn,se);
    e2.CapSize=0;
    e2.LineWidth=2;
    e2.Color = rgb('grey');
    
end



