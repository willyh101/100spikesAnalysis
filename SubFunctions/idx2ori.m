function oris = idx2ori(stimIdxs, orientations)
% can include nans, just has to be the same length

assert(numel(unique(stimIdxs)) == numel(orientations), ...
    'Number of unique stimIdxs does not equal number of orientations.')

oris = arrayfun(@(x) orientations(x), stimIdxs);