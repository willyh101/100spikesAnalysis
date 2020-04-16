function [All, outVars] = redCellTuningAnalysis(All, outVars, opts)

redTuningPlot = [];
redOSIplot = [];
notRedTuningPlot = [];
notRedOSIplot = [];
isVisRedPlot = [];
isVisOtherPlot = [];
redCurvesPlot = [];
otherCurvesPlot = [];

alpha = opts.visAlpha;

for ind = 1:numel(All)
    
    oris = [nan 0:45:315];
    isRed = All(ind).out.red.isRed;
    isVisRed = All(ind).out.anal.pVisR < alpha & isRed;
    isVisRedCells = find(isVisRed);
    isVisOther = All(ind).out.anal.pVisR < alpha & ~isRed;
    isVisOtherCells = find(isVisOther);
    
    All(ind).out.red.isVis = isVisRed;
    All(ind).out.red.isVisCells = isVisRedCells;
    
    disp(['Found ' num2str(numel(isVisRedCells)) ' visually responsive and red cells.'])
    
    % red cells first
    % tuning idx
    redTuningIdx = All(ind).out.anal.prefOri(isVisRed);
    All(ind).out.red.redTuningIdx = redTuningIdx;    
    % real ori
    redTuningOri = idx2ori(redTuningIdx, [nan 0:45:315]); % actual oris
    All(ind).out.red.redTuningOri = redTuningOri;
    % osi
    redOSI = All(ind).out.anal.osi(isVisRed);
    All(ind).out.red.redOSI = redOSI;
    % curves
    redCurves = All(ind).out.anal.oriCurve(:, isVisRed);
    All(ind).out.red.redCurves = redCurves;
    redCurvesSEM = All(ind).out.anal.oriCurveSEM(:, isVisRed);
    All(ind).out.red.redCurvesSEM = redCurvesSEM;
    
    % other cells
    % pref cond
    notRedTuningIdx = All(ind).out.anal.prefOri(~isVisOther);
    All(ind).out.red.notRedTuningIdx = notRedTuningIdx;
    % ori
    notRedTuningOri = idx2ori(notRedTuningIdx, [nan 0:45:315]);
    All(ind).out.red.notRedTuningOri = notRedTuningOri;
    % osi
    notRedOSI = All(ind).out.anal.osi(~isVisOther);
    All(ind).out.red.notRedOSI = notRedOSI;
    % curves
    otherCurves = All(ind).out.anal.oriCurve(:, isVisOther);
    All(ind).out.red.otherCurves = otherCurves;
    otherCurvesSEM = All(ind).out.anal.oriCurveSEM(:, isVisOther);
    All(ind).out.red.otherCurvesSEM = otherCurvesSEM;
    
    % for plotting...
    redTuningPlot = [redTuningPlot redTuningOri];
    redOSIplot = [redOSIplot redOSI];
    notRedTuningPlot = [notRedTuningPlot notRedTuningOri];
    notRedOSIplot = [notRedOSIplot notRedOSI];
    isVisRedPlot = [isVisRedPlot isVisRed];
    isVisOtherPlot = [isVisOtherPlot isVisOther];
    
    % save into outVars
    outVars.redCellTuning{ind} = redTuningOri;
    outVars.redOSI{ind} = redOSI;
    outVars.notRedTuning{ind} = notRedTuningOri;
    outVars.notRedOSI{ind} = notRedOSI;
    outVars.redCurves{ind} = redCurves;
    outVars.redCurvesSEM{ind} = redCurvesSEM;
    outVars.otherCurves{ind} = otherCurves;
    outVars.otherCurvesSEM{ind} = otherCurvesSEM;
    
end

%% comparison of pref ori and OSI for red and not red
figure(40)
clf

subplot(2,2,1)
histogram(redTuningPlot, 9, 'FaceColor', 'r')
title('Red Cell Tuning')
ylabel('Count')
xlabel('Preferred Orientation')

subplot(2,2,2)
histogram(redOSIplot, 20, 'FaceColor', 'r')
title('Red Cell Selectivity')
ylabel('Count')
xlabel('OSI')

subplot(2,2,3)
histogram(notRedTuningPlot, 9, 'FaceColor', rgb('gray'))
title('Not Red Cell Tuning')
ylabel('Count')
xlabel('Preferred Orientation')

subplot(2,2,4)
histogram(notRedOSIplot, 20, 'FaceColor', rgb('gray'))
title('Not Red Cell Selectivity')
ylabel('Count')
xlabel('OSI')

%% %%
% [All, outVars] = getTuningCurve(All, opts, outVars);

%% example cells for each
h41 = figure(41);
clf

num_example = opts.numExamples;

% clear redcurves redcurvesSEM notredcurves notredcurvesSEM
% for ind = 1:numel(All)
% 
%     redMask = find(All(ind).out.red.isRed & All(ind).out.anal.osi > osi_thresh & All(ind).out.anal.pVisR < alpha);
%     notRedMask = find(~All(ind).out.red.isRed & All(ind).out.anal.osi > osi_thresh & All(ind).out.anal.pVisR < alpha);
% 
%     redcurves{ind} = All(ind).out.anal.oriCurve(:, redMask);
%     redcurvesSEM{ind} = All(ind).out.anal.oriCurveSEM(:, redMask);
% 
%     notredcurves{ind} = All(ind).out.anal.oriCurve(:, notRedMask);
%     notredcurvesSEM{ind} = All(ind).out.anal.oriCurveSEM(:, notRedMask);
% 
% end

% redcurves = redcurves{:};
% redcurvesSEM = redcurvesSEM{:};
% notredcurves = notredcurves{:};
% notredcurvesSEM = notredcurvesSEM{:};
% %


red2plt = horzcat(outVars.redCurves{:});
red2pltSEM = horzcat(outVars.redCurvesSEM{:});
other2plt = horzcat(outVars.otherCurves{:});
other2pltSEM = horzcat(outVars.otherCurvesSEM{:});

try
    assert(num_example <= size(red2plt,2))
catch
    num_example = size(red2plt,2);
    warning(['Number of requested examples (red cells) is greater than number of cells available. Setting to ' num2str(num_example) '.'])
end


redex = randperm(size(red2plt, 2), num_example);
red2plt = red2plt(:, redex);
red2pltSEM = red2pltSEM(:, redex);

otherex = randperm(size(other2plt, 2), num_example);
other2plt = other2plt(:, otherex);
other2pltSEM = other2pltSEM(:, otherex);

nplts = num_example;

c = 0;
for i = 1:nplts
    c = c + 1;
    subplot(2, nplts, c)
    
    data = red2plt(:, i);
    err = red2pltSEM(:, i);
    
    e = errorbar(oris, data, err);
    e.Color = 'red';
    e.LineWidth = 1;
    ylabel('zdf')
    xlabel('ori')
end

for i = 1:nplts
    c = c + 1;
    subplot(2, nplts, c)
    
    data = other2plt(:, i);
    err = other2pltSEM(:, i);
    
    e = errorbar(oris, data, err);
    e.Color = rgb('grey');
    e.LineWidth = 1;
    ylabel('zdf')
    xlabel('ori')
end