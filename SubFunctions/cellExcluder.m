function [All, results] = cellExcluder(All,opts); 

if ~isfield(opts, 'peakThreshold')
    opts.peakThreshold =0;
end

numExp = numel(All);

for ind = 1:numExp
    dat = All(ind).out.exp.dataToUse;
    sz=size(dat);
    dat = reshape(dat,[sz(1) sz(2)*sz(3)]);
    meanDat = mean(dat,2);
    maxDat = max(dat,[],2); 
    
    cellsToExclude = meanDat<opts.minMeanThreshold | meanDat>opts.maxMeanThreshold...
        | maxDat<opts.peakThreshold;
    results{ind}= cellsToExclude;
    if opts.verbose
        disp(['ind: ' num2str(ind) '. ' num2str(sum(cellsToExclude)) ' Excluded, ' num2str(mean(cellsToExclude)*100,2) '%'])
    end
    
    All(ind).out.anal.cellsToExclude = cellsToExclude';
    All(ind).out.anal.cellsToInclude = ~cellsToExclude';

end