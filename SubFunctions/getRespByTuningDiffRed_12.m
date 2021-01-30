function outVars = getRespByTuningDiffRed_12(All, outVars)

ensemblesToUseList = find(outVars.ensemblesToUse );

OSIList = outVars.meanEnsOSI(ensemblesToUseList);
[sortedOSI, sortOrder] = sort(OSIList);
ensemblesToUseList=ensemblesToUseList(sortOrder); 

diffRange = -15:30:330; %0:33.75:315;

% osiToUse = outVars.meanEnsOSI;
oriToUse = outVars.meanEnsOri;

figure(137);clf;
mRespByOriDiff=[];
for i=1:numel(ensemblesToUseList)
    
    
    ens = ensemblesToUseList(i);
    ind = outVars.ensIndNumber(ens);
    hNum = outVars.ensHNumber(ens);
    
    if ~isfield(All(ind).out, 'red')
        continue
    end
    
    thisOSI = outVars.meanEnsOSI(ens);
    thisOri = oriToUse(ens);
    
    cellToUse = ~All(ind).out.anal.offTargetRisk(hNum-1,:) & ...
        ~All(ind).out.anal.ROIinArtifact' & ...
        All(ind).out.anal.pVisR < 0.05 & ...
        All(ind).out.red.isRed;
    
    oris = outVars.prefOris{ind};
    oriOptions  = [NaN 0:30:330];
    oris = oriOptions(oris); 
    
   
    oriDiff = abs(thisOri-oris);
    oriDiff(oriDiff>180) = abs((thisOri-360)-(oris(oriDiff>180)));
    oriDiff(oriDiff>180) = abs((thisOri)-(oris(oriDiff>180)-360));
% oriDiff = (thisOri-oris);
% oriDiff(oriDiff>180) = ((thisOri-360)-(oris(oriDiff>180)));
% oriDiff(oriDiff>180) = ((thisOri)-(oris(oriDiff>180)-360));
% oriDiff(oriDiff<-180) = ((thisOri-360)-(oris(oriDiff<-180)));
% oriDiff(oriDiff<-180) = ((thisOri)-(oris(oriDiff<-180)-360));
    
    respToPlot = outVars.mRespEns{ens}(cellToUse);
    oriDiffToPlot = oriDiff(cellToUse);
    
    plot(oriDiffToPlot,respToPlot,'.')
    title(['Ens: ' num2str(ens) ' OSI: ' num2str(thisOSI)]);
    xlim([-5 185]);
    
    tempOriDiff =nan([1,numel(diffRange)-1]);
    for k=1:numel(diffRange)-1
        tempOriDiff(k) = nanmean(respToPlot(oriDiffToPlot>=diffRange(k) & oriDiffToPlot<diffRange(k+1)));
    end
    hold on
    p=plot(diffRange(1:end-1),tempOriDiff);
    p.Marker = 'o';
    p.LineWidth =2;
    p.Color = 'r';
    
%     pause(0.1)

    
    mRespByOriDiffRed(i,:) = tempOriDiff;
    
end

outVars.sortedOSIred = sortedOSI;
outVars.mRespByOriDiffRed = mRespByOriDiffRed;