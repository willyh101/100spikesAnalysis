%% Example cells for distance

ensHNumber = outVars.ensHNumber;
ensIndNumber = outVars.ensIndNumber;
distBins = opts.distBins;
ensemblesToUse = outVars.ensemblesToUse;

ensemblesToUse = ensemblesToUse &...
                 outVars.ensOSI > 0.67 &...
                 outVars.numCellsEachEns == 10;

disp(['num ensembles: ' num2str(sum(ensemblesToUse))])
             
% binRange = [0 50];
binRange = [0 50];

% distToUse = distBins(binToUse);

% choose or randomly select ensembles
% ensToUse = randsample(find(ensemblesToUse),1);
ensToUse = 641;

hNum = ensHNumber(ensToUse); % ensembleID
ind = ensIndNumber(ensToUse); % All out ID

ensPO = outVars.ensPO(ensToUse);

stimDistance =  All(ind).out.anal.StimDistance;
ROIinArtifact = All(ind).out.anal.ROIinArtifact;
offTargetRisk = All(ind).out.anal.offTargetRisk;
prefs = idx2ori(All(ind).out.anal.prefOri, [nan 0:45:315]);
pVisR = All(ind).out.anal.pVisR;
tgs = All(ind).out.anal.targets;
targeted = zeros(size(pVisR));
targeted(tgs) = 1;

h = All(ind).out.exp.stimParams.roi{hNum};
tg = All(ind).out.exp.rois{h};
dists = stimDistance(tg,:);

minStimDistance = min(dists,[],1);



cellsToUse = ~ROIinArtifact' ...
             & ~targeted ...
             & ~offTargetRisk(h,:) ...
             & minStimDistance > binRange(1) ...
             & minStimDistance <= binRange(2) ...
             & pVisR < 0.05 ...
             ...& prefs == ensPO ...
             ;

disp(['num cells:   ' num2str(sum(cellsToUse))])
cs = find(cellsToUse);
%%this plots the average TS of all cells for that condition

boot = 0; % bootstrap the number of cells?

if boot > 0
    c = randsample(cs,ncellToPlot);
else
    c = cs;
end
v=1; % no vis stim cond

figure(311)
clf

us = unique(All(ind).out.exp.stimID);
s = us(hNum);

trialsToUse = All(ind).out.exp.lowMotionTrials &...
              All(ind).out.exp.stimSuccessTrial &...
              ismember(All(ind).out.exp.stimID, s) &...
              ismember(All(ind).out.exp.visID, 1);

dat = squeeze(nanmean(All(ind).out.exp.zdfData(c,:,trialsToUse),3));
dat = dat-nanmean(dat(:, 1:5),2);



% fillPlot is data,[],lineCol,edgeCol,faceCol,fAlpha
[hL2, hF2] = fillPlot(dat,[], 'sem', rgb('black'), [], rgb('royalblue'), .5);
title('Cell Average')

tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
xlabel('Time')
ylabel('ZDF')
suptitle(['Ensemble ' num2str(ensToUse)])



    
%% this plots 4 individual cells

% choose or randomly select cells from eligible
cs_sample = randsample(cs,4);
v=1; % no vis stim cond

figure(310)
clf

% % do the sort
% bwin = 1:5;
% mwin = 7:12;
% 
% mdat = squeeze(mean(mean(All(ind).out.exp.zdfData(cellsToUse,mwin,trialsToUse),2),3));
% mbdat = squeeze(mean(mean(All(ind).out.exp.zdfData(cellsToUse,bwin,trialsToUse),2),3));
% mdat = mdat - mbdat;

for i=1:4
    subplot(2,2,i)
    c = cs_sample(i);
    us = unique(All(ind).out.exp.stimID);
    s = us(hNum);
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
                  All(ind).out.exp.stimSuccessTrial &...
                  ismember(All(ind).out.exp.stimID, s) &...
                  ismember(All(ind).out.exp.visID, 1);
    
    dat = squeeze(All(ind).out.exp.zdfData(c,:,trialsToUse));
    dat = dat-mean(dat(1:5,:),1);
    
    % fillPlot is data,[],lineCol,edgeCol,faceCol,fAlpha
    [hL2, hF2] = fillPlot(dat',[], 'sem', rgb('black'), [], rgb('royalblue'), .5);
    title(['Cell ' num2str(c) ', ' num2str(round(minStimDistance(c))), ' um'])
    
    tcks = 0:6:30;
    xticks(tcks)
    xticklabels((tcks/6)-1)
    xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
    xlabel('Time')
    ylabel('ZDF')
    
end
    
suptitle(['Ensemble ' num2str(ensToUse)])

%% plot a specific pair of cells

figure(312)
clf
% 
% ensToUse = 624;
% c1 = 410;
% c2 = 371;
% ensToUse = 640;
ensToUse = 643;
c1 = 577;
c2 = 53;


hNum = ensHNumber(ensToUse); % ensembleID
ind = ensIndNumber(ensToUse); % All out ID

stimDistance =  All(ind).out.anal.StimDistance;
ROIinArtifact = All(ind).out.anal.ROIinArtifact;
offTargetRisk = All(ind).out.anal.offTargetRisk;
pVisR = All(ind).out.anal.pVisR;
tgs = All(ind).out.anal.targets;
targeted = zeros(size(pVisR));
targeted(tgs) = 1;

h = All(ind).out.exp.stimParams.roi{hNum};
tg = All(ind).out.exp.rois{h};
dists = stimDistance(tg,:);

minStimDistance = min(dists,[],1);


subplot(1,2,1)
c=c1;
us = unique(All(ind).out.exp.stimID);
s = us(hNum);
trialsToUse = All(ind).out.exp.lowMotionTrials &...
              All(ind).out.exp.stimSuccessTrial &...
              ismember(All(ind).out.exp.stimID, s) &...
              ismember(All(ind).out.exp.visID, 1);
          
dat = squeeze(All(ind).out.exp.zdfData(c,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);

[hL2, hF2] = fillPlot(dat',[], 'sem', rgb('black'), [], rgb('grey'), .5);

hold on

ax = gca();
ax.FontSize = 12;
ax.LineWidth = 1;
x1 = ax.YLim(1);
x2 = abs(x1) + ax.YLim(2);



r = rectangle('position',[6 x1 6 x2]);
r.FaceColor = [rgb('FireBrick') 0.25];
r.LineStyle = 'none';
uistack(r,'bottom')

yline(0,'--k', 'LineWidth',1.5)

% text(9,ax.YLim(2)*0.9,{'Ensemble','Stim'},'HorizontalAlignment','center')

% title(['Close - ' num2str(round(minStimDistance(c))) ' \mum'])
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xlabel('Time (s)')
ylabel('\Delta z-scored dF/F')
box off

c=c2;
subplot(1,2,2)
us = unique(All(ind).out.exp.stimID);
s = us(hNum);
trialsToUse = All(ind).out.exp.lowMotionTrials &...
              All(ind).out.exp.stimSuccessTrial &...
              ismember(All(ind).out.exp.stimID, s) &...
              ismember(All(ind).out.exp.visID, 1);
          
dat = squeeze(All(ind).out.exp.zdfData(c,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);


[hL2, hF2] = fillPlot(dat',[], 'sem', rgb('black'), [], rgb('grey'), .5);

hold on

ax = gca();
ax.FontSize = 12;
ax.LineWidth = 1;
x1 = ax.YLim(1);
x2 = abs(x1) + ax.YLim(2);


r = rectangle('position',[6 x1 6 x2]);
r.FaceColor = [rgb('FireBrick') 0.25];
r.LineStyle = 'none';
uistack(r,'bottom')

yline(0,'--k', 'LineWidth',1.5)

% text(9,ax.YLim(2)*0.9,{'Ensemble','Stim'},'HorizontalAlignment','center')

% title(['Far - ' num2str(round(minStimDistance(c))) ' \mum'])
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xlabel('Time (s)')
ylabel('\Delta z-scored dF/F')
box off

%% Plot 4 cells to 2 different stimuli

figure(313)
clf
% 
% ensToUse = 624;
% c1 = 410;
% c2 = 371;
% ensToUse = 640;
e1 = 640;
e2 = 643;
c1 = 577;
c2 = 53;

ensToUse = e1;
hNum = ensHNumber(ensToUse); % ensembleID
ind = ensIndNumber(ensToUse); % All out ID

stimDistance =  All(ind).out.anal.StimDistance;
ROIinArtifact = All(ind).out.anal.ROIinArtifact;
offTargetRisk = All(ind).out.anal.offTargetRisk;
pVisR = All(ind).out.anal.pVisR;
tgs = All(ind).out.anal.targets;
targeted = zeros(size(pVisR));
targeted(tgs) = 1;

h = All(ind).out.exp.stimParams.roi{hNum};
tg = All(ind).out.exp.rois{h};
dists = stimDistance(tg,:);

minStimDistance = min(dists,[],1);


subplot(2,2,1)
c=c1;
us = unique(All(ind).out.exp.stimID);
s = us(hNum);
trialsToUse = All(ind).out.exp.lowMotionTrials &...
              All(ind).out.exp.stimSuccessTrial &...
              ismember(All(ind).out.exp.stimID, s) &...
              ismember(All(ind).out.exp.visID, 1);
          
dat = squeeze(All(ind).out.exp.zdfData(c,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);

[hL2, hF2] = fillPlot(dat',[], 'sem', rgb('black'), [], rgb('grey'), .5);

hold on

ax = gca();
ax.FontSize = 12;
ax.LineWidth = 1;
x1 = ax.YLim(1);
x2 = abs(x1) + ax.YLim(2);



r = rectangle('position',[6 x1 6 x2]);
r.FaceColor = [rgb('FireBrick') 0.25];
r.LineStyle = 'none';
uistack(r,'bottom')

yline(0,'--k', 'LineWidth',1.5)

% text(9,ax.YLim(2)*0.9,{'Ensemble','Stim'},'HorizontalAlignment','center')

% title(['Close - ' num2str(round(minStimDistance(c))) ' \mum'])
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xlabel('Time (s)')
ylabel('\Delta z-scored dF/F')
box off

c=c2;
subplot(2,2,2)
us = unique(All(ind).out.exp.stimID);
s = us(hNum);
trialsToUse = All(ind).out.exp.lowMotionTrials &...
              All(ind).out.exp.stimSuccessTrial &...
              ismember(All(ind).out.exp.stimID, s) &...
              ismember(All(ind).out.exp.visID, 1);
          
dat = squeeze(All(ind).out.exp.zdfData(c,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);


[hL2, hF2] = fillPlot(dat',[], 'sem', rgb('black'), [], rgb('grey'), .5);

hold on

ax = gca();
ax.FontSize = 12;
ax.LineWidth = 1;
x1 = ax.YLim(1);
x2 = abs(x1) + ax.YLim(2);


r = rectangle('position',[6 x1 6 x2]);
r.FaceColor = [rgb('FireBrick') 0.25];
r.LineStyle = 'none';
uistack(r,'bottom')

yline(0,'--k', 'LineWidth',1.5)

% text(9,ax.YLim(2)*0.9,{'Ensemble','Stim'},'HorizontalAlignment','center')

% title(['Far - ' num2str(round(minStimDistance(c))) ' \mum'])
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xlabel('Time (s)')
ylabel('\Delta z-scored dF/F')
box off


%%%---SECOND ENSEMBLE---%%%

ensToUse = e2;
hNum = ensHNumber(ensToUse); % ensembleID
ind = ensIndNumber(ensToUse); % All out ID

stimDistance =  All(ind).out.anal.StimDistance;
ROIinArtifact = All(ind).out.anal.ROIinArtifact;
offTargetRisk = All(ind).out.anal.offTargetRisk;
pVisR = All(ind).out.anal.pVisR;
tgs = All(ind).out.anal.targets;
targeted = zeros(size(pVisR));
targeted(tgs) = 1;

h = All(ind).out.exp.stimParams.roi{hNum};
tg = All(ind).out.exp.rois{h};
dists = stimDistance(tg,:);

minStimDistance = min(dists,[],1);


subplot(2,2,3)
c=c1;
us = unique(All(ind).out.exp.stimID);
s = us(hNum);
trialsToUse = All(ind).out.exp.lowMotionTrials &...
              All(ind).out.exp.stimSuccessTrial &...
              ismember(All(ind).out.exp.stimID, s) &...
              ismember(All(ind).out.exp.visID, 1);
          
dat = squeeze(All(ind).out.exp.zdfData(c,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);

[hL2, hF2] = fillPlot(dat',[], 'sem', rgb('black'), [], rgb('grey'), .5);

hold on

ax = gca();
ax.FontSize = 12;
ax.LineWidth = 1;
x1 = ax.YLim(1);
x2 = abs(x1) + ax.YLim(2);



r = rectangle('position',[6 x1 6 x2]);
r.FaceColor = [rgb('FireBrick') 0.25];
r.LineStyle = 'none';
uistack(r,'bottom')

yline(0,'--k', 'LineWidth',1.5)

% text(9,ax.YLim(2)*0.9,{'Ensemble','Stim'},'HorizontalAlignment','center')

% title(['Close - ' num2str(round(minStimDistance(c))) ' \mum'])
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xlabel('Time (s)')
ylabel('\Delta z-scored dF/F')
box off

c=c2;
subplot(2,2,4)
us = unique(All(ind).out.exp.stimID);
s = us(hNum);
trialsToUse = All(ind).out.exp.lowMotionTrials &...
              All(ind).out.exp.stimSuccessTrial &...
              ismember(All(ind).out.exp.stimID, s) &...
              ismember(All(ind).out.exp.visID, 1);
          
dat = squeeze(All(ind).out.exp.zdfData(c,:,trialsToUse));
dat = dat-mean(dat(1:5,:),1);


[hL2, hF2] = fillPlot(dat',[], 'sem', rgb('black'), [], rgb('grey'), .5);

hold on

ax = gca();
ax.FontSize = 12;
ax.LineWidth = 1;
x1 = ax.YLim(1);
x2 = abs(x1) + ax.YLim(2);


r = rectangle('position',[6 x1 6 x2]);
r.FaceColor = [rgb('FireBrick') 0.25];
r.LineStyle = 'none';
uistack(r,'bottom')

yline(0,'--k', 'LineWidth',1.5)

% text(9,ax.YLim(2)*0.9,{'Ensemble','Stim'},'HorizontalAlignment','center')

% title(['Far - ' num2str(round(minStimDistance(c))) ' \mum'])
tcks = 0:6:30;
xticks(tcks)
xticklabels((tcks/6)-1)
xlabel('Time (s)')
ylabel('\Delta z-scored dF/F')
box off




