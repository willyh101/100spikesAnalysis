%% Example cells for distance

ensHNumber = outVars.ensHNumber;
ensIndNumber = outVars.ensIndNumber;
distBins = opts.distBins;

ensemblesToUse = ensemblesToUse &...
                 outVars.ensOSI > 0.67 &...
                 outVars.numCellsEachEns == 10;

sum(ensemblesToUse)
             
binRange = [10 30];
% binRange = [75 125];

% distToUse = distBins(binToUse);

% choose or randomly select ensembles
% ensToUse = randsample(find(ensemblesToUse),1);
ensToUse = 624;

hNum = ensHNumber(ensToUse); % ensembleID
ind = ensIndNumber(ensToUse); % All out ID

stimDistance =  All(ind).out.anal.StimDistance;
ROIinArtifact = All(ind).out.anal.ROIinArtifact;
offTargetRisk = All(ind).out.anal.offTargetRisk;
pVisR = All(ind).out.anal.pVisR;

h = All(ind).out.exp.stimParams.roi{hNum};
tg = All(ind).out.exp.rois{h};
dists = stimDistance(tg,:);

minStimDistance = min(dists,[],1);



cellsToUse = ~ROIinArtifact' ...
             & ~offTargetRisk(h,:) ...
             & minStimDistance > binRange(1) ...
             & minStimDistance <= binRange(2) ...
             & pVisR < 0.05 ...
             ;

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

dat = squeeze(mean(All(ind).out.exp.zdfData(c,:,trialsToUse),3));
dat = dat-mean(dat(:, 1:5),2);

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

