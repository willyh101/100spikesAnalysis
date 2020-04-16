function oris = idx2ori(stimIdxs, orientations)
% can include nans, just has to be the same length

oris = arrayfun(@(x) orientations(x), stimIdxs);