function [All outVars] = defineCorrelationTypesOnVis(All, outVars)
numExps = numel(All);

for ind = 1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
            numCells = All(ind).out.anal.numCells; 

    if ~isfield(All(ind).out, 'vis')
        fprintf(' No Vis... ')
        
        
        All(ind).out.anal.SpontCorr = nan([numCells numCells]);
    All(ind).out.anal.SpCoP = nan([numCells numCells]);
    
    All(ind).out.anal.AllCorr = nan([numCells numCells]);
    All(ind).out.anal.AlCoP = nan([numCells numCells]);
    
    All(ind).out.anal.AllMCorr = nan([numCells numCells]);
    All(ind).out.anal.AmCoP =nan([numCells numCells]);
    
    All(ind).out.anal.SignalCorr = nan([numCells numCells]);
    All(ind).out.anal.SiCoP =nan([numCells numCells]);
    
    All(ind).out.anal.NoiseCorr = nan([numCells numCells]);
    All(ind).out.anal.NoCoP = nan([numCells numCells]);
    else
    %Spont Corr - correlation coefficient on time series from no stim
    %period
    trialsToUse = All(ind).out.vis.lowMotionTrials &...
        All(ind).out.vis.lowRunTrials &...
        All(ind).out.vis.visID == 1;
    unrollData = All(ind).out.vis.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    if isempty(unrollData)
        SpontCorr  = nan([numCells numCells]);
        SpCoP = nan([numCells numCells]);
    else
        [SpontCorr SpCoP] = corr(unrollData');
    end
    
    %AllCorr - the correlation coef on all time series
    trialsToUse = All(ind).out.vis.lowMotionTrials &...
        All(ind).out.vis.lowRunTrials ;
      
    unrollData = All(ind).out.vis.zdfData(:,:,trialsToUse);
    sz = size(unrollData);
    unrollData = reshape(unrollData,[sz(1) sz(2)*sz(3)]);
    
    [AllCorr AlCoP] = corr(unrollData');
    
    %All corr mean - correlation coef of response (not time series)
    trialsToUse = All(ind).out.vis.lowMotionTrials &...
        All(ind).out.vis.lowRunTrials;
       
    unrollData = All(ind).out.vis.rdData(:,trialsToUse);
    sz = size(unrollData);
    
    [AllMCorr AmCoP] = corr(unrollData');
    
    %noise corr - correlation coef of residual trial response (not time
    %series) i.e. trial response - mean response for that condition
    trialsToUse = All(ind).out.vis.lowMotionTrials &...
        All(ind).out.vis.lowRunTrials;
       
    vs = unique(All(ind).out.vis.visID);
    vs(vs==0)=[];
    
    unrollData = [];
    meanResps = [];
    for k = 1:numel(vs)
        v = vs(k);
        trialsToUseThis = trialsToUse & All(ind).out.vis.visID==v;
        
        dataPart = All(ind).out.vis.rdData(:,trialsToUseThis);
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
    end
    
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
