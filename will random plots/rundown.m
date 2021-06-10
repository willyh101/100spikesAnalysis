figure(77)
clf

out = All(1).out;
% dfsource = 'zdfData';
dfsource = 'allData';
% 
% if strcmp(dfsource, 'zdfData')
%     dflbl = 'zdf';
% else
%     dflbl = 'df/f';
%     [dfData, ~] = computeDFFwithMovingBaseline(out.vis.allData);
%     out.vis.dfData = dfData;
%     [dfData, ~] = computeDFFwithMovingBaseline(out.exp.allData);
%     out.exp.dfData = dfData;
% end

stim = 2;
vis = [1 2]; % these should be the blank conditions

us = unique(out.exp.stimID);
s = us(stim);

vs = unique(out.exp.visID);
v = vs(vis);

h = out.exp.stimParams.roi{stim};
tg = out.exp.holoTargets{h};
tg(isnan(tg))=[];


%%--- ChroME rundown figure ---%%

trialsToUse = ismember(out.exp.stimID, s) &...
              ismember(out.exp.visID, v) &...
              out.exp.lowMotionTrials &...
              out.exp.stimSuccessTrial;

dat = squeeze(out.exp.(dfsource)(tg,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);

ntrials = sum(trialsToUse);

resp = mean(dat(12:18,:),1);

scatter(find(trialsToUse), resp, 'filled', 'MarkerFaceColor','blue')

hold on
ft = polyfit(find(trialsToUse), resp, 1);
ft_xs = 1:length(trialsToUse);
ftvals = polyval(ft, ft_xs);
plot(ft_xs, ftvals, 'LineWidth', 1, 'Color', 'blue')

 
%%----- 2nd stim -----%%
stim = 3;
vis = [1 2]; % these should be the blank conditions

us = unique(out.exp.stimID);
s = us(stim);

vs = unique(out.exp.visID);
v = vs(vis);

h = out.exp.stimParams.roi{stim};
tg = out.exp.holoTargets{h};
tg(isnan(tg))=[];




trialsToUse = ismember(out.exp.stimID, s) &...
              ismember(out.exp.visID, v) &...
              out.exp.lowMotionTrials &...
              out.exp.stimSuccessTrial;

dat = squeeze(out.exp.(dfsource)(tg,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);

ntrials = sum(trialsToUse);

resp = mean(dat(12:18,:),1);

scatter(find(trialsToUse), resp, 'filled', 'MarkerFaceColor','red')

ft = polyfit(find(trialsToUse), resp, 1);
ft_xs = 1:length(trialsToUse);
ftvals = polyval(ft, ft_xs);
plot(ft_xs, ftvals, 'LineWidth', 1, 'Color', 'red')

xlabel('Trial Number')
ylabel('z-scored df/f')
title('ChroME rundown?')

figure(78)

subplot(1,2,1)

