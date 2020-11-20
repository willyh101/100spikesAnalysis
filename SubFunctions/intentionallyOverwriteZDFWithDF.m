function All = intentionallyOverwriteZDFWithDF(All)

warning off backtrace

warning(['You are replacing the ZDF data with DF data in the All.out ' ...
        'struct! However, his will not alter the contents of the outfiles '...
        'unless of course you re-save them. Proceed with caution. Also, '...
        'any experiments lacking dfData will be removed from the All(out) '...
        'struct entirely, thus breaking any indexing into the loadList, '...
        'so be sure to run this AFTER allLoadListErrorFixer.'])
    
numExpts = numel(All);
c=0;

for ind = 1:numExpts
    try
        % overwrite it
        All(ind).out.exp.zdfData = All(ind).out.exp.dfData;
    catch ME
        % if there isn't a dfData field, skip that expt and remove it
        if strcmp(ME.identifier,'MATLAB:nonExistentField')
            c = c+1;
            disp(['No dfData for experiment ' num2str(ind) '. Skipping and removing from All(out).'])
            exptsNoDF(c) = ind; % to be removed
        else
            % if we don't know what happened we want to raise the error
            rethrow(ME)
        end  
    end
end
try
    All(exptsNoDF) = []; % removed from All(out)
catch
    disp('All expts had dfData. That''s great.')
end

disp(['There are now ' num2str(numel(All)) ' experiments.'])