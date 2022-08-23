function stateN(ensemblesToUse,outVars)
ensemblesToUse = logical(ensemblesToUse); 
ensIndNumber = outVars.ensIndNumber;
names = outVars.names; 

disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);