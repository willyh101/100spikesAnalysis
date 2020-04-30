

ensemblesToUse = outVars.ensemblesToUse;

popRespAll = outVars.allEnsResp(ensemblesToUse)';
popRespIso = outVars.isoEnsResp(ensemblesToUse)';
popRespOrtho = outVars.orthoEnsResp(ensemblesToUse)';
popRespNotCoTuned = outVars.notCoTunedResp(ensemblesToUse)';
popRespNotVis = outVars.notVisEnsResp(ensemblesToUse)';

ncells = outVars.numCellsEachEns(ensemblesToUse)';
osi = outVars.ensOSI(ensemblesToUse)';
size = outVars.ensMeaD(ensemblesToUse)';
corr = outVars.ensAlCo(ensemblesToUse)';

names = {'popAll', 'popIso', 'popOrtho', 'popNotCoTuned', 'popNotVis', ...
           'numCells', 'osi', 'size', 'corr'};

ensTable = table(popRespAll, popRespIso, popRespOrtho, popRespNotCoTuned, ...
            popRespNotVis, ncells, osi, size, corr, 'VariableNames', names);
        

grps = {'popAll', 'popIso', 'popOrtho', 'popNotCoTuned', 'popNotVis'};     
means = grpstats(ensTable, 'numCells', 'mean', grps);

figure(420)
hold on
for g = grps
    p = plot(m.numCells, means{:, {strcat('mean_', g{:})}});
    p.LineWidth = 1;
end
    

