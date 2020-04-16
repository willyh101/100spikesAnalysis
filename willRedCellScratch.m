%% Red cell co-tuning section

redTuningPlot = [];
redOSIplot = [];
notRedTuningPlot = [];
notRedOSIplot = [];

for ind = 1:numel(All)
    oris = [nan 0:45:315];
    isRed = All(ind).out.red.isRed;
    
    % red cells first
    redTuningIdx = All(ind).out.anal.prefOri(isRed);
    All(ind).out.red.redTuningIdx = redTuningIdx;

    redTuningOri = arrayfun(@(x) oris(x), redTuningIdx);
    All(ind).out.red.redTuningOri = redTuningOri;

    redOSI = All(ind).out.anal.osi(isRed);
    All(ind).out.red.redOSI = redOSI;
    
    % other cells
    notRedTuningIdx = All(ind).out.anal.prefOri(~isRed);
    All(ind).out.red.notRedTuningIdx = notRedTuningIdx;
    
    notRedTuningOri = arrayfun(@(x) oris(x), notRedTuningIdx);
    All(ind).out.red.notRedTuningOri = notRedTuningOri;

    notRedOSI = All(ind).out.anal.osi(~isRed);
    All(ind).out.red.notRedOSI = notRedOSI;
    
    % for plotting...
    redTuningPlot = [redTuningPlot redTuningOri];
    redOSIplot = [redOSIplot redOSI];
    notRedTuningPlot = [notRedTuningPlot notRedTuningOri];
    notRedOSIplot = [notRedOSIplot notRedOSI];
end

% comparison of pref ori and OSI for red and not red
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
%%
[All, outVars] = getTuningCurve(All, opts, outVars);
% [All, outVars] = calcOSI(All, outVars);
%%
% example cells for each
figure(41)
clf

num_example = 5;
osi_thresh = 0.3;

% redex = datasample(find( ...
%     All(2).out.red.isRed & ...
%     All(2).out.anal.osi > osi_thresh & ...
%     All(2).out.anal.pVisR < 0.05), num_example, 'Replace', false);

redex = datasample(find( ...
    All(2).out.red.isRed & ...
    All(2).out.anal.osi > osi_thresh), num_example, 'Replace', false);

notredex = datasample(find( ...
    ~All(2).out.red.isRed & ...
    All(2).out.anal.osi > osi_thresh & ...
    All(2).out.anal.pVisR < 0.05), num_example, 'Replace', false);

curves = All(2).out.anal.oriCurve;
curvesSEM = All(2).out.anal.oriCurveSEM;

redcurves = curves(:, redex);
redcurvesSEM = curvesSEM(:, redex);

notredcurves = curves(:, notredex);
notredcurvesSEM = curvesSEM(:, notredex);


nplts = size(redcurves, 2);

c = 0;
for i = 1:nplts
    c = c + 1;
    subplot(2, nplts, c)
    
    data = redcurves(:, i);
    err = redcurvesSEM(:, i);
    
    e = errorbar(oris, data, err);
    e.Color = 'red';
    e.LineWidth = 1;
end

for i = 1:nplts
    c = c + 1;
    subplot(2, nplts, c)
    
    data = notredcurves(:, i);
    err = notredcurvesSEM(:, i);
    
    e = errorbar(oris, data, err);
    e.Color = rgb('grey');
    e.LineWidth = 1;
end

gcf;
axis tight
    
    












