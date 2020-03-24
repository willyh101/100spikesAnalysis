%% table demo

dataIn = All(1).out;
visTable = makeVisTable(dataIn);
head(visTable) % show the top of table

%% indexing/selecting

a = visTable(visTable.visID==1,:);
b = visTable(:,{'visID', 'lowRun'});
c = visTable(visTable{:,'visID'}==1,:);

trialsToUse = visTable((visTable.lowRun==true & visTable.lowMotion==true),:);

%% grouping

% summaries, calculate means and stuff in table form, usually the best way
groupsummary(trialsToUse, 'visID', 'mean')
groupsummary(trialsToUse, 'visID', {'mean', 'max'}, 'zdf')
means = groupsummary(trialsToUse, 'visID', {'mean','min'}, 'zdf');
vals = means.mean_zdf; % extracts the mean values out of the table

% on multiple groupings
groupsummary(visTable,{'visID','lowRun'},{'mean','median'}, 'zdf')

% add back to the original table
joinedback = join(trialsToUse, means);

%% sub-select
trialsToUse = visTable(visTable.lowMotion==1 & visTable.lowRun==1,:);
head(trialsToUse)

