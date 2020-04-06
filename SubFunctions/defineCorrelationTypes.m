function [All outVars] = defineCorrelationTypes(All, outVars)
numExps = numel(All);

for ind = 1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    us = unique(All(ind).out.exp.stimID);
    
    %Spont Corr - correlation coefficient on time series from no stim
    %period
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) &...
        All(ind).out.exp.visID == 1;
    unrollData = All(ind).out.exp.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    [SpontCorr SpCoP] = corr(unrollData');
    
    %AllCorr - the correlation coef on all time series
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) ;
    unrollData = All(ind).out.exp.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    [AllCorr AlCoP] = corr(unrollData');
    
    %All corr mean - correlation coef of response (not time series)
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) ;
    unrollData = All(ind).out.exp.rdData(:,trialsToUse);
    sz = size(unrollData);
    
    [AllMCorr AmCoP] = corr(unrollData');
    
    %noise corr - correlation coef of residual trial response (not time
    %series) i.e. trial response - mean response for that condition
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial &...
        All(ind).out.exp.stimID == us(1) ;
    vs = unique(All(ind).out.exp.visID);
    unrollData = [];
    meanResps = [];
    for k = 1:numel(vs)
        v = vs(k);
        trialsToUseThis = trialsToUse & All(ind).out.exp.visID==v;
        
        dataPart = All(ind).out.exp.rdData(:,trialsToUseThis);
        mData = mean(dataPart');
        meanResps(k,:) =  mData;
        dataPart = dataPart-mData';
        
        unrollData = cat(2,unrollData,dataPart);
    end
    
    [NoiseCorr NoCoP] = corr(unrollData');
    [SignalCorr SiCoP] = corr(meanResps);
    
    
    All(ind).out.anal.SpontCorr = SpontCorr;
    All(ind).out.anal.SpCoP = SpCoP;
    
    All(ind).out.anal.AllCorr = AllCorr;
    All(ind).out.anal.AlCoP = AlCoP;
    
    All(ind).out.anal.AllMCorr = AllMCorr;
    All(ind).out.anal.AmCoP =AmCoP;
    
    All(ind).out.anal.SignalCorr = SignalCorr;
    All(ind).out.anal.SiCoP =SiCoP;
    
    All(ind).out.anal.NoiseCorr = NoiseCorr;
    All(ind).out.anal.NoCoP = NoCoP;
    
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% Determine the Ensemble CoCorrelation

ensSpCo=[];ensAlCo=[];ensAmCo=[];ensSiCo=[];ensNoCo=[];
for ind = 1:numExps
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    
    clear ensembleSpCo ensembleAlCo ensembleAmCo ensembleSiCo ensembleNoCo
    %     for i =1:numel(All(ind).out.exp.holoTargets)
    %         ht = All(ind).out.exp.holoTargets{i};
    %         ht(isnan(ht))=[];
    c=0;
    for i= 1:numStims
        holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
        if holo>0
            c=c+1;
            ht = All(ind).out.exp.holoTargets{holo};
            ht(isnan(ht))=[];
            
            corrToUse = All(ind).out.anal.SpontCorr;
            corMat = corrToUse(ht,ht);
            corMat(logical(eye(numel(ht))))=nan;
            ensembleSpCo(c) = nanmean(corMat(:));
            
            corrToUse = All(ind).out.anal.AllCorr;
            corMat = corrToUse(ht,ht);
            corMat(logical(eye(numel(ht))))=nan;
            ensembleAlCo(c) =nanmean(corMat(:));
            
            corrToUse = All(ind).out.anal.AllMCorr;
            corMat = corrToUse(ht,ht);
            corMat(logical(eye(numel(ht))))=nan;
            ensembleAmCo(c) = nanmean(corMat(:));
            
            corrToUse = All(ind).out.anal.SignalCorr;
            corMat = corrToUse(ht,ht);
            corMat(logical(eye(numel(ht))))=nan;
            ensembleSiCo(c) =nanmean(corMat(:));
            
            corrToUse = All(ind).out.anal.NoiseCorr;
            corMat = corrToUse(ht,ht);
            corMat(logical(eye(numel(ht))))=nan;
            ensembleNoCo(c) =nanmean(corMat(:));
        end
    end
    All(ind).out.anal.ensembleSpCo = ensembleSpCo;
    All(ind).out.anal.ensembleAlCo = ensembleAlCo;
    All(ind).out.anal.ensembleAmCo = ensembleAmCo;
    All(ind).out.anal.ensembleSiCo = ensembleSiCo;
    All(ind).out.anal.ensembleNoCo = ensembleNoCo;
    
    ensSpCo = cat(2,ensSpCo,ensembleSpCo);
    ensAlCo = cat(2,ensAlCo,ensembleAlCo);
    ensAmCo = cat(2,ensAmCo,ensembleAmCo);
    ensSiCo = cat(2,ensSiCo,ensembleSiCo);
    ensNoCo = cat(2,ensNoCo,ensembleNoCo);
end

outVars.ensSpCo=ensSpCo;
outVars.ensAlCo=ensAlCo;
outVars.ensAmCo=ensAmCo;
outVars.ensSiCo=ensSiCo;
outVars.ensNoCo=ensNoCo;

ensemblesToUse = outVars.ensemblesToUse;
figure(214);clf
dat = {ensSpCo(ensemblesToUse), ensAlCo(ensemblesToUse), ensAmCo(ensemblesToUse),  ensSiCo(ensemblesToUse), ensNoCo(ensemblesToUse)};
names = {'Spont' 'All' 'All (v2)' 'Signal' 'Noise'};
fancyPlotSpread(dat,names);
title('Ensemble Mean Correlations by type')
ylabel('Correlation (Rho)')

%%plot Pop Response by Correlation
f3 = figure(203);
clf(f3);

popResponseEns = outVars.popResponseEns; 
numCellsEachEns = outVars.numCellsEachEns;

for i=1:5
    subplot(5,1,i)
    dataToUse = dat{i};
    scatter(dataToUse,popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
    % scatter(1:sum(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
    
    title([names{i} ' Correlation'])
    
    xlabel(['Correlation of Ensemble'])
    ylabel('Population Mean Response')
    % title('OSIs by Ensemble Size')
    set(gcf(),'Name','OSIs by Ensemble Size')
    cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
    cb.Label.String = 'Number of Cells in Ensemble';
    r = refline(0);
    r.LineStyle =':';
end
