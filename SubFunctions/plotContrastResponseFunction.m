function [outVars] = plotContrastResponseFunction(All,outVars,opts)
%%
%% Contrast response functions
distBinCRF      = opts. distBinCRF ;
visAlphaCRF     = opts.visAlphaCRF;

ensIndNumber    = outVars.ensIndNumber;
visPercent      = outVars.visPercent;
ensemblesToUse  = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns; 

ensSizes = unique(numCellsEachEns(ensemblesToUse)); 
colorList = colorMapPicker(numel(ensSizes),outVars.defaultColorMap); 

numExps = numel(All); 

counter = 0;
clear CRF CRFNS CRFDist CRFNSDist
for ind = 1:numExps
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    pVisR = All(ind).out.anal.pVisR;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    us = unique(All(ind).out.exp.stimID);
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    for i=1:numel(us)
        s = us(i);
        h = All(ind).out.exp.stimParams.roi{i};
        
        if h==0
            continue
        end
        cellsToUse = ~ROIinArtifact' & pVisR<visAlphaCRF &...
            ~any(offTargetRisk);%~offTargetRisk(h,:); %~any(offTargetRisk);% (h,:); 
  
        
        counter=counter+1;
        CRF(counter,:) = nanmean(respMat(i,:,cellsToUse) - baseMat(i,:,cellsToUse),3);
        CRFNS(counter,:) = nanmean(respMat(1,:,cellsToUse) - baseMat(1,:,cellsToUse),3);
        
        for d=1:numel(distBinCRF)-1
            cellsToUseDist = cellsToUse &...
                All(ind).out.anal.minDistbyHolo(h+1,:) <=distBinCRF(d+1) &...
                All(ind).out.anal.minDistbyHolo(h+1,:) >distBinCRF(d) ;
            
            CRFDist(counter,:,d) = nanmean(respMat(i,:,cellsToUseDist) - baseMat(i,:,cellsToUseDist),3);
            CRFNSDist(counter,:,d) = nanmean(respMat(1,:,cellsToUseDist) - baseMat(1,:,cellsToUseDist),3);
        
        end
    end
end

%catch missing data
shouldBeNAN = CRF==0 & CRFNS ==0;
CRF(shouldBeNAN) = nan;
CRFNS(shouldBeNAN) = nan;       
CRFDiff = CRF-CRFNS;

shouldBeNAN = CRFDist==0 & CRFNSDist ==0;
CRFDist(shouldBeNAN) = nan;
CRFNSDist(shouldBeNAN) = nan;       
CRFDiffDist = CRFDist-CRFNSDist;

highVisPercentIndMoreStringent = ~ismember(ensIndNumber,find(visPercent<0.1)); %remove low vis responsive experiments

ensemblesToUseThis = ensemblesToUse;% & highVisPercentIndMoreStringent ;%& ensIndNumber==8;% highVisPercentIndMoreStringent ;

figure(13);clf
subplot(1,2,1)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat1 = CRF(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:);
meanDat = nanmean(dat1);
stdDat = nanstd(dat1);
numpDat = sum(~isnan(dat1));
semDat = stdDat./sqrt(numpDat);

data1{i}=dat1;

hold on
errorbar(meanDat,semDat,'linewidth',2,'color',colorList{i})

% dat2 = CRFNS(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:);

dat2 = CRFNS(ensemblesToUseThis ,:);
indsOrigin{i} = ensIndNumber(ensemblesToUseThis & numCellsEachEns==ensSizes(i));
data2{i}=dat2;
meanDat = nanmean(dat2);
stdDat = nanstd(dat2);
numpDat = sum(~isnan(dat2));
semDat = stdDat./sqrt(numpDat);
e2=errorbar(meanDat,semDat,'linewidth',2,'color',rgb('grey'),'linestyle',':');
% e2.Color = colorList{i};

end

ylabel('Z-score dF')
xticks(1:6);
xticklabels({'0%' '1%' '4%' '10%' '40%' '100%'})
xlabel('Contrast')

subplot(1,2,2)

for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
dat1 = CRFDiff(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:);
meanDat = nanmean(dat1);
stdDat = nanstd(dat1);
numpDat = sum(~isnan(dat1));
semDat = stdDat./sqrt(numpDat);

data1{i}=dat1;

hold on
errorbar(meanDat,semDat,'linewidth',2,'color',colorList{i})
end
r = refline(0);
r.LineWidth = 2;
r.LineStyle =':';
r.Color = rgb('grey');

ylabel('\Delta Z-Scored dF')
xlabel('Contrast')
xticks(1:6);
xticklabels({'0%' '1%' '4%' '10%' '40%' '100%'})

clear ax
figure(14);clf
for c=1:numel(distBinCRF)-1
    subplot(1,numel(distBinCRF)-1,c)
for i = 1:size(ensSizes,2)
% subplot(1,size(ensSizes),i)
% dat1 = CRFDiffDist(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:,c);
dat1 = CRFDist(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:,c);

meanDat = nanmean(dat1);
stdDat = nanstd(dat1);
numpDat = sum(~isnan(dat1));
semDat = stdDat./sqrt(numpDat);

data1{i}=dat1;

hold on
errorbar(meanDat,semDat,'linewidth',2,'color',colorList{i})

% dat2 = CRFNSDist(ensemblesToUseThis & numCellsEachEns==ensSizes(i),:,c);

dat2 = CRFNSDist(ensemblesToUseThis ,:,c);
indsOrigin{i} = ensIndNumber(ensemblesToUseThis & numCellsEachEns==ensSizes(i));
data2{i}=dat2;
meanDat = nanmean(dat2);
stdDat = nanstd(dat2);
numpDat = sum(~isnan(dat2));
semDat = stdDat./sqrt(numpDat);
% errorbar(meanDat,semDat,'linewidth',2,'color',rgb('grey'),'linestyle',':')

end


ylabel('Z-score dF')
xticks(1:6);
xticklabels({'0%' '1%' '4%' '10%' '40%' '100%'})
xlabel('Contrast')

title([num2str(distBinCRF(c)) ' to ' num2str(distBinCRF(c+1)) ' \mum from Target'])
end
linkaxes