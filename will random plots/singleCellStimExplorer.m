% out = All(2).out;
% dfsource = 'dfData';
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

stim = 6;
vis = [1]; % these should be the blank conditions

us = unique(out.exp.stimID);
s = us(stim);

vs = unique(out.exp.visID);
v = vs(vis);

h = out.exp.stimParams.roi{stim};
tg = out.exp.holoTargets{h};
tg(isnan(tg))=[];




%%------ stimmed cell figure --------%%
% first, plot tuning curve, holo response, and vis resonse for the stimmed cell 
figure(66)
clf
sgtitle(['Stim ID ' num2str(stim) ', ROI ' num2str(tg)])

% holo stim response
subplot(1,3,1)

trialsToUse = ismember(out.exp.stimID, s) &...
              ismember(out.exp.visID, v) &...
              out.exp.lowMotionTrials &...
              out.exp.stimSuccessTrial;

dat = squeeze(out.exp.(dfsource)(tg,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);
ntrials = sum(trialsToUse);

[hL, hF] = fillPlot(dat',[],'ci');
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
xlabel('Time')
ylabel(dflbl)
title('Holo stim response (40p, 40 Hz)')

% vis response curve
subplot(1,3,2)

trialsToUse = out.vis.visID > 1 &...
              out.vis.lowMotionTrials;
          
dat = squeeze(out.vis.(dfsource)(tg,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);
ntrials = size(dat,2);

[hL, hF] = fillPlot(dat',[],'ci');
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
xlabel('Time')
ylabel(dflbl)
title('Mean vis stim response, ori epoch')

% tuning curve
subplot(1,3,3)
e = errorbar(out.anal.oriCurve(2:end,tg), out.anal.oriCurveSEM(2:end,tg));
e.LineWidth = 2;
e.Color = 'k';
title(['Ori Tuning Curve, OSI=' num2str(out.anal.osi(tg))])
ylabel(dflbl)
xticks(1:12)
xticklabels(0:30:330)
xtickangle(45)
xlabel('Direction')
po = out.anal.prefOri(tg)-1;
lbl_y = e.YData(po) + 1.1*e.YPositiveDelta(po);
text(po, lbl_y, '*', 'Color','r','FontSize',30)




%%------- other cells figure --------%%
f67 = figure(67);
clf
rws = 1;
cls = 4;

sgtitle(['Stim ID ' num2str(stim) ', ROI ' num2str(tg)])



% mean pop response all others
subplot(rws,cls,1)

trialsToUse = ismember(out.exp.stimID, s) &...
              ismember(out.exp.visID, v) &...
              out.exp.lowMotionTrials &...
              out.exp.stimSuccessTrial;
          
cellsToUse = ~out.anal.ROIinArtifact' &...
             ~out.anal.offTargetRisk(stim-1,:) &...
             out.anal.pVisR < 0.05;

isred = out.red.isRed(cellsToUse);

disp(['Total blue cells:     ' num2str(sum(out.red.isRed))])
disp(['Blue cells vis resp:  ' num2str(sum(out.red.isRed & out.anal.pVisR < 0.05))])
disp(['Blue cells used:      ' num2str(sum(isred))])
         
dat = squeeze(mean(out.exp.(dfsource)(cellsToUse,:,trialsToUse),3));
dat = dat-mean(dat(:,1:5),2);
% popdat = squeeze(mean(dat,1));
% ntrials = size(dat,2);

% fillPlot is data,[],lineCol,edgeCol,faceCol,fAlpha

[hL2, hF2] = fillPlot(dat(isred,:),[], 'ci', rgb('red'), [], rgb('red'), 0.5);
hold on
[hL1, hF1] = fillPlot(dat(~isred,:),[],'ci', rgb('black'), [], rgb('black'), 0.5);

tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
xlabel('Time')
ylabel(dflbl)
title('Mean Population response')
legend([hL1,hL2], {'Pyramidal Cells','SST Cells'})



% scatter plot of all cells - no stim condition
subplot(rws,cls,2)

trialsToUse = out.exp.stimID == us(1) &...
              ismember(out.exp.visID, [1 2]) &...
              out.exp.lowMotionTrials;
          
dat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
dat = dat - bdat;

% dat = sort(dat);
s1 = scatter(find(~isred),dat(~isred), 'filled');
s1.MarkerFaceAlpha = 0.7;
s1.MarkerFaceColor = rgb('grey');
hold on
s2 = scatter(find(isred),dat(isred), 'filled');
s2.MarkerFaceColor = rgb('lightsalmon');

yline(0, '--','LineWidth',1)
ylabel(dflbl)
% ax = gca();
% ax.XAxis.Visible = 'off';
xlabel('Cells')
title({'Individual cell responses','No stim'})
legend([s1,s2], {'Pyramidal Cells','SST Cells'})
ax(1) = gca();


% scatter plot of all cells
subplot(rws,cls,3)

trialsToUse = ismember(out.exp.stimID, s) &...
              ismember(out.exp.visID, v) &...
              out.exp.lowMotionTrials &...
              out.exp.stimSuccessTrial;
          
dat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
dat = dat - bdat;

% dat = sort(dat);
s1 = scatter(find(~isred),dat(~isred), 'filled');
s1.MarkerFaceAlpha = 0.7;
% s1.MarkerFaceColor = rgb('grey');
hold on
s2 = scatter(find(isred),dat(isred), 'filled');
s2.MarkerFaceColor = 'r';

yline(0, '--','LineWidth',1)
ylabel(dflbl)
% ax = gca();
% ax.XAxis.Visible = 'off';
xlabel('Cells')
title('Individual cell responses')
legend([s1,s2], {'Pyramidal Cells','SST Cells'})
ax(2) = gca();
linkaxes(ax)



% % indiv cells by tuning (errorbar)
% subplot(rws,cls,3)
% colormap(f67, 'viridis')
% oris = 0:30:330;
% pref = oris(po);
% all_pref = arrayfun(@(x) oris(x), out.anal.prefOri(out.anal.prefOri>1)-1);
% diffToPO = mod(all_pref - pref,360);
% dat2 = mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2);
% error
% 
% 
% s1 = scatter(1:sum(isred),sort(dat(isred)),out.red.redOSI*50, diffToPO, 'filled');
% c = colorbar;
% c.Label.String = 'PO difference';
% c.Ticks = 0:30:330;
% ylabel(dflbl)
% xlabel('SST cells (sorted)')
% title('SST response by difference to stim PO')




% indiv cells by tuning (scatter since there aren't many)
subplot(rws,cls,4)
colormap(f67, 'viridis')
oris = 0:30:330;
pref = oris(po);
diffToPO = mod(out.red.redTuningOri - pref,360);
s1 = scatter(1:sum(isred),sort(dat(isred)),out.red.redOSI*100, diffToPO, 'filled');
c = colorbar;
c.Label.String = 'PO difference';
c.Ticks = 0:30:330;
ylabel(dflbl)
xlabel('SST cells (sorted)')
title('SST response by difference to stim PO')













