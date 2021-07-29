    8                    \
    
    
    
    
    
    l;,.'$wrnknumEns = numel(outVars.mRespEns);

oriVals = [NaN 0:30:330];

clear temp temp2

for i=1:numEns
    diffsPossible = [0 30 60 90 120 150 180];
    indNum = outVars.ensIndNumber(i);
    cellOris = oriVals(outVars.prefOris{indNum});
    cellOrisDiff = abs(cellOris-outVars.ensPO(i));
    cellOrisDiff(cellOrisDiff>180) = abs(cellOrisDiff(cellOrisDiff>180)-360);
    temp{i} = cellOrisDiff;
    
    hNum = outVars.ensHNumber(i);
    
    cellsToUse = ~All(indNum).out.anal.ROIinArtifact' &...
                 ~All(indNum).out.anal.offTargetRisk(hNum-1,:) &...
                 All(indNum).out.anal.pVisR < 0.05;
             
    temp2{i} = cellsToUse';
end

outVars.cellOrisDiff = temp;
outVars.cellsToUse = temp2;
    

%%
ensData = outVars.mRespEns(ensemblesToUse);
ensDists = outVars.distToEnsemble(ensemblesToUse);
oriDiff = outVars.cellOrisDiff(ensemblesToUse);
cellsToUse = outVars.cellsToUse(ensemblesToUse);
    
    
ensData = vertcat(ensData{:});
ensDists = horzcat(ensDists{:});
oriDiff = horzcat(oriDiff{:});
cellsToUse = vertcat(cellsToUse{:});

ensData = ensData(cellsToUse);
ensDists = ensDists(cellsToUse);
oriDiff = oriDiff(cellsToUse);

%% DISTANCE

[sorted_dists, sort_idx] = sort(ensDists);

sorted_resps = ensData(sort_idx);

figure(760)
clf

s1 = scatter(sorted_dists, sorted_resps', 'filled');
s1.MarkerFaceAlpha = 0.5;
xlabel('Distance')
ylabel('z-scored \DeltaF/F')

%% TUNING
figure(780)
clf

clear diffr_resp diffr_sem
for i=1:numel(diffsPossible)
    
    diffr = diffsPossible(i);
    
    diffr_resp(i) = nanmean(ensData(oriDiff==diffr));
    diffr_sem(i) = sem(ensData(oriDiff==diffr));
end

e1 = errorbar(diffr_resp, diffr_sem);
e1.LineWidth = 2;
ylabel('zscored \DeltaF/F')
xlabel('\Delta PO')
xticklabels(0:30:180)

    
    
    