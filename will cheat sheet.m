%%-----Will's Handy MATLAB snippets-----%%

% go from condition -> true value (stimID -> condition, visID -> condition)
% for each value in idx_arry, use that as an index into another arr (vals)
idx_arr = [1 2 4 4 3 2];
vals = ['cond1', 'cond2', 'cond3', 'cond4'];
result = arrayfun(@(x) vals(x), idx_arr);

% group an array by another variable
% ex: grouping ensemble response by mouse, arrays must be same length
grouping = [1 1 1 2 2];
input_arr = [0.13, 0.11, 0.2, -0.4, -1.1];
unique_groups = unique(grouping);
n_per_cond = diff([0, find(diff(unique_groups)), numel(unique_groups)]);
result = cellfun(@cell2mat, mat2cell(input_arr, 1, n_per_cond), 'un', 0);

% how to use ismember...
result = ismember(long_arr, arr_of_conds);

% select a field from a struct with a string
field = 'hologram';
result = some_structure.(field);
