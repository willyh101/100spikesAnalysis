function plotCompareRedCellScatter(outVars)

%%---RED CELLS---%%

ensemblesToUse = outVars.ensemblesToUse;
popRespCoTuned = outVars.coTunedRedEnsResp(ensemblesToUse);
popRespOrthoTuned = outVars.orthoTunedRedEnsResp(ensemblesToUse);
popRespNotCoTuned = outVars.notCoTunedRedEnsResp(ensemblesToUse);
popRespNotVisResp = outVars.notVisRespRedEnsResp(ensemblesToUse);

order = 1:numel(popRespCoTuned);
eOSI = outVars.ensOSI(ensemblesToUse);
dist = outVars.ensMeaD(ensemblesToUse);
allcorr = outVars.ensAlCo(ensemblesToUse);
sigcorr = outVars.ensSiCo(ensemblesToUse);
size = outVars.numCellsEachEns(ensemblesToUse);


compareBy = [order; eOSI; dist; allcorr; sigcorr; size];

compareNames = {'Order of Being Done' 'Ensemble OSI' 'Ensemble Mean Distance' 'Ensemble All Corr' 'Ensemble Signal Corr' 'Ensemble Size'};
    
fig = figure(333);
clf
colormap(fig, 'viridis')
hold on

for ind=1:numel(compareNames)
    subplot(2,3,ind)

    x = compareBy(ind,:);
    xname = compareNames(ind);
    
    s1 = scatter(x, popRespCoTuned, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on
    s2 = scatter(x, popRespOrthoTuned, 'filled', 'MarkerFaceAlpha', 0.7);
    s3 = scatter(x, popRespNotCoTuned, 'filled', 'MarkerFaceAlpha', 0.7);
    s4 = scatter(x, popRespNotVisResp, 'filled', 'MarkerFaceAlpha', 0.7);
    
    title({'Mean Population Response', 'PV Cells'})
    xlabel(xname)
    ylabel('PV Cell Population Mean Response')
    r = refline(0);
    r.LineStyle = ':';
    
    if ind==1
        legend([s1, s2, s3, s4], {'Co-Tuned', 'Ortho-Tuned', 'Not Co-Tuned', sprintf('Not Vis\nResponsive')})
    end
end

%%---OTHER CELLS---%%

ensemblesToUse = outVars.ensemblesToUse;
popRespCoTuned = outVars.coTunedOtherEnsResp(ensemblesToUse);
popRespOrthoTuned = outVars.orthoTunedOtherEnsResp(ensemblesToUse);
popRespNotCoTuned = outVars.notCoTunedOtherEnsResp(ensemblesToUse);
popRespNotVisResp = outVars.notVisRespOtherEnsResp(ensemblesToUse);


fig = figure(334);
clf
colormap(fig, 'viridis')
hold on

for ind=1:numel(compareNames)
    subplot(2,3,ind)

    x = compareBy(ind,:);
    xname = compareNames(ind);
    
    s1 = scatter(x, popRespCoTuned, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on
    s2 = scatter(x, popRespOrthoTuned, 'filled', 'MarkerFaceAlpha', 0.7);
    s3 = scatter(x, popRespNotCoTuned, 'filled', 'MarkerFaceAlpha', 0.7);
    s4 = scatter(x, popRespNotVisResp, 'filled', 'MarkerFaceAlpha', 0.7);
    
    title({'Mean Population Response','Pyr Cells'})
    xlabel(xname)
    ylabel('Pyramidal Cell Population Mean Response')
    r = refline(0);
    r.LineStyle = ':';
    
    if ind==1
        legend([s1, s2, s3, s4], {'Co-Tuned', 'Ortho-Tuned', 'Not Co-Tuned', sprintf('Not Vis\nResponsive')})
    end
end
    
    
    
    
    
    