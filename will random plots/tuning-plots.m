%% example tuning curves vis responses etc
figure(2)
clf

sgtitle('Example Tuning Curve - TRE-ChroME2s')

nplots = 16;
vis_cells = find(outVars.isVisR{1});
r = randperm(numel(vis_cells), nplots);
cells = vis_cells(r);
oris = [0:45:315];

for i=1:nplots
    subplot(4,4,i)
    c = cells(i);
    curve = outVars.tuningCurves{1};
    curve = curve(2:9,c);
    curveSEM = outVars.tuningCurvesSEM{1};
    curveSEM = curveSEM(2:9,c);
    e = errorbar(oris, curve, curveSEM);
    e.LineWidth = 1;
    title(['Cell  ' num2str(c)])
    xticks(oris)
    xtickangle(45)
end
%% OSI and PO plots

osis = outVars.osi{1};
osis = osis(outVars.isVisR{1});
osis = osis(~isnan(osis));

pos = outVars.prefOris{1};
pos = pos(outVars.isVisR{1});
pos = pos(pos>1);

figure(3)
clf

subplot(1,2,1)
histogram(osis,10)
xlabel('OSI')
ylabel('Count')

subplot(1,2,2)
histogram(pos,8)
xticks(2:10)
xticklabels(oris)
xtickangle(45)
xlabel('Preferred Direction')
ylabel('Count')