% generic pattern to extract data from response matrices for a single
% outfile

% this assumes you are dropping in a outfile, not loading them into the
% All(out).ind structure

% this is just ver1, in order to get the tuning of ensembles, you would need
% to load up the files into the All(out) data structure to do some initial
% analysis... I will follow-up with ver2 that can do all of that. so for
% now this only works on the vis data


% first, choose your data source
% you can change to out.vis.dfData if it exists
dataSource = out.vis.zdfData;


% trials to use section...

% change visID into orientations
% nan is grey screen
vs = out.vis.visID;
orisUsed = [nan 0:45:315];
oris = arrayfun(@(x) orisUsed(x), vs);

% select which oris to look at here
oriToUse = [0 180];

trialsToUse = ismember(oris, oriToUse);


% cells to use section

% get stimmed cells from exp section
cellsToUse = out.exp.targetedCells;
cellsToUse = cellsToUse(~isnan(cellsToUse))';


% result is cells x time x trials
data = dataSource(cellsToUse, :, trialsToUse);