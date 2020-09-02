function plotResponseByDifferenceinAnglePref(outVars,All, opts)

ensemblesToUseList = find(outVars.ensemblesToUse );

OSIList = outVars.meanEnsOSI(ensemblesToUseList);
[sortedOSI, sortOrder] = sort(OSIList);
ensemblesToUseList=ensemblesToUseList(sortOrder); 

diffRange = -22.5:45:315; %0:33.75:315;

% osiToUse = outVars.meanEnsOSI;
oriToUse = outVars.meanEnsOri;

figure(137);clf;
mRespByOriDiff=[];
for i=1:numel(ensemblesToUseList)
    if mod(i,50)==0
        disp('.')
    else
        fprintf('.')
    end
    clf
    ens = ensemblesToUseList(i);
    ind = outVars.ensIndNumber(ens);
    hNum = outVars.ensHNumber(ens);
    
    thisOSI = outVars.meanEnsOSI(ens);
    thisOri = oriToUse(ens);
    
    cellToUse = ~All(ind).out.anal.offTargetRisk(hNum-1,:) & ...
        ~All(ind).out.anal.ROIinArtifact' & ...
        All(ind).out.anal.pVisR < 0.05 ;
    
    oris = outVars.prefOris{ind};
    oriOptions  = [NaN 0:45:315];
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

    
    mRespByOriDiff(i,:) = tempOriDiff;
    
end
%%
figure(138);clf
% imagesc(mRespByOriDiff);
% colormap rdbu
% caxis([-0.1 0.1])
subplot(1,2,1)
meanByOriDIff = (nanmean(mRespByOriDiff));
semByOriDiff = nanstd(mRespByOriDiff)./sqrt(sum(~isnan(mRespByOriDiff)));
e1 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDIff,semByOriDiff);
e1.LineWidth = 2;

% plot only good OSIs
hold on
goodOSIthreshold = opts.goodOSIthresh;
datToPlot = mRespByOriDiff(sortedOSI>=goodOSIthreshold,:);
meanByOriDIff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
e2 = errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDIff,semByOriDiff);
e2.LineWidth = 2;



ylim([-0.05 0.01]) 
xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')
xticks(0:45:135)
legend({'All OSIs', ['OSI > ' num2str(goodOSIthreshold)]})


subplot(1,2,2)
% imToPlot = mRespByOriDiff;
% imToPlot(isnan(imToPlot))=0;
% imagesc(imToPlot)
% caxis([-0.2 0.2])
% xlim([0.5 7.5])
% pT1 = prctile(sortedOSI,25);
% pT2 = prctile(sortedOSI,50);
% pT3 = prctile(sortedOSI,75);

numEns = numel(sortedOSI); 
hold on
% datToPlot = mRespByOriDiff(1:round(numEns/4),:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(round(numEns/4):round(numEns/2),:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(round(numEns/2):round(numEns*3/4),:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)
% 
% datToPlot = mRespByOriDiff(round(numEns*3/4):end,:);
% meanByOriDIff = (nanmean(datToPlot));
% semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
% errorbar(diffRange(1:end-1),meanByOriDIff,semByOriDiff)


datToPlot = mRespByOriDiff(sortedOSI<0.25,:);
meanByOriDIff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDIff,semByOriDiff)

datToPlot = mRespByOriDiff(sortedOSI>=0.25 & sortedOSI<=0.75,:);
meanByOriDIff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDIff,semByOriDiff)

datToPlot = mRespByOriDiff(sortedOSI>0.75,:);
meanByOriDIff = (nanmean(datToPlot));
semByOriDiff = nanstd(datToPlot)./sqrt(sum(~isnan(datToPlot)));
errorbar(diffRange(1:end-1)-diffRange(1),meanByOriDIff,semByOriDiff)


% ylim([-0.05 0.01])
 xlim([-5 185])
xlabel('\Delta Preferred Angle (Deg)')
ylabel('\Delta Response')

legend('Low OSI <0.25','Mid OSI','High OSI >0.75')



