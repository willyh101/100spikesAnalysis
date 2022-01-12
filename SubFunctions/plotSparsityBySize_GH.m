function plotSparsityBySize_GH(All,outVars)

subtractBaseline =1;

ensemblesToUseList = find(outVars.ensemblesToUse );
for i=1:numel(ensemblesToUseList)
    if mod(i,50)==0
        disp('.')
    else
        fprintf('.')
    end
    
    ens = ensemblesToUseList(i);
    ind = outVars.ensIndNumber(ens);
    hNum = outVars.ensHNumber(ens);
    
    try
        cellToUse = ~All(ind).out.anal.offTargetRisk(hNum-1,:) ...
        & ~All(ind).out.anal.ROIinArtifact' ...
         & All(ind).out.anal.pVisR < 0.01 ...
        & ~All(ind).out.red.isRed...
        ;
    catch
    cellToUse = ~All(ind).out.anal.offTargetRisk(hNum-1,:)...
        & ~All(ind).out.anal.ROIinArtifact' ...
        &  All(ind).out.anal.pVisR < 0.01 ...
        ;
    end
    
    
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
end

ensSize = outVars.numCellsEachEns(ensemblesToUseList);

uEns = unique(ensSize);

for i=1:numel(uEns);
    data{i}= allSparse(ensSize==uEns(i));
    names{i}= string(uEns(i));
end

figure(55);clf
cmap = colorMapPicker(numel(uEns),outVars.defaultColorMap);
% cmap{end+1} = rgb('grey');

p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

ylabel('Sparsity L2/L1')
xlabel('Ensemble Size')