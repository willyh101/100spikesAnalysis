%% loadlists

ppsfLoadList

% loadPath = 'T:\Outfiles';
loadPath = 'E:\outfiles';
addpath(genpath(loadPath))

%% Load data

numExps = numel(loadList);
disp(['There are ' num2str(numExps) ' Exps in this LoadList'])
if numExps ~= 0
    clear All
    if ~iscell(loadList)
        numExps=1;
        temp = loadList;
        clear loadList;
        loadList{1} = temp;
    end

    for ind = 1:numExps
        pTime =tic;
        fprintf(['Loading Experiment ' num2str(ind) '...']);
        All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
        fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
    end
else
    disp('Did you press this by accident?')
end

%%error fixer
%CAUTION! ERRORS WILL OCCUR IF YOU RUN MORE THAN ONCE!
[All] = allLoadListErrorFixer(All,loadList);

%%set Data To use
for ind=1:numExps
    All(ind).out.exp.dataToUse = All(ind).out.exp.dfData;
end
disp('Data To Use is set')

% %% clean Data, and create fields.
% 
opts.FRDefault=6;
opts.recWinRange = [0.5 1.5]; %[0.5 1.5];[1.5 2.5];%[0.5 1.5];% %from vis Start in s [1.25 2.5];


%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull

%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded
opts.runValPercent = 0.75; %percent of frames that need to be below run threshold

[All, outVars] = cleanData(All,opts);

%% rematch to the holos we care about only

stims = [1 25];
depthList = [0 30 55];
clear targs
for ind=1:2
    Mapping = zeros([20 3]);
    cc = 0;
    for s=stims
%         r = All(ind).out.exp.rois{s};
        r = All(ind).out.exp.holoRequest.rois{s};
        sCoM = All(ind).out.exp.holoRequest.targets(r,1:2);
        sDepth = uint16(All(ind).out.exp.holoRequest.targets(r,3));
        for i=1:numel(sDepth)
            sDepth(i) = find(sDepth(i)==depthList);
        end
        sCoM = sCoM + All(ind).out.exp.offsets;

        for j=1:size(sCoM,1)
            d = sDepth(j);
            Mapping(j,1) = d;
            a = All(ind).out.exp.allCoM;
            a(All(ind).out.exp.allDepth~=d,:) = nan;

            cc = cc+1;
            b = sCoM(j,:);
            [Mapping(cc,3), Mapping(cc,2)] = min(sqrt(sum((a-b).^2,2)));
        end
    end
    tg_cells = Mapping(:,2);
    tg_cells(Mapping(:,3)>10)=nan;
    targs{ind} = tg_cells;
end

%% calculate responses

origHoloIdx = [1:24; 25:48];
cellIdx = [1:10; 11:20];

c=0;

clear resp
for ind=1:2
    for h=1:2
        c=c+1;
        hIdx = origHoloIdx(h,1);
        htg = targs{ind};
        r = htg(cellIdx(h,:));
        r = r(~isnan(r));
        cellsToUse = r;

        all_us = All(ind).out.exp.uniqueStims;
        stim_range = origHoloIdx(h,:);
        us = all_us(stim_range+1);
        
        clear temp
        for i=1:numel(us)
            s = us(i);

            trialsToUse = All(ind).out.exp.stimID == s & ...
                All(ind).out.exp.lowMotionTrials;
            hnum = All(ind).out.exp.stimParams.roi{i+hIdx};
            holo_shift = All(ind).out.exp.rois{hnum};
            
            % compute responses
            dat = All(ind).out.exp.dataToUse;
            win = All(ind).out.anal.recWinUsed;
            bwin = All(ind).out.anal.bwinToUse;

            cellVal = mean(mean(dat(cellsToUse,win(1):win(2),trialsToUse),2),3)';
            bVal = mean(mean(dat(cellsToUse,bwin(1):bwin(2),trialsToUse),2),3)';
            cellVal = cellVal-bVal;
            temp(i,:) = cellVal;
%             err_temp(i,:) =
            % error
            mcellVal = squeeze(mean(dat(cellsToUse,win(1):win(2),trialsToUse),2));
            mbVal = squeeze(mean(dat(cellsToUse,bwin(1):bwin(2),trialsToUse),2));
            result = mcellVal - mbVal;
            mVal = mean(result,2);
%             eVal = 

        end
        resp{c} = temp;
    end
end

resp = double(cell2mat(resp)');
% ax and rad are cell x cond
ax = resp(:,1:12);
rad = resp(:,13:24);

%% calculate hwhms

stim_success = 0.25;


x = [-3:3:30];
z = [-6:6:60];


clear hwr hwa fits
for i=1:size(rad,1)
    % radial
    y = rad(i,:);
    try
        [ff, gof] = fit(x',y', 'gauss1');
        hwr(i) = hwhm(ff.c1);
        fits(i).radial.fit = ff;
        fits(i).radial.gof = gof;
        fits(i).radial.yvals = y;
        fits(i).radial.xvals = x;
        fits(i).radial.xvals_adj = x-ff.b1;
        fits(i).radial.hwhm = hwr(i);
    catch
        hwr(i) = nan;
        fits(i).radial = nan;
    end

    % axial
    y = ax(i,:);
    try
        [ff, gof] = fit(z',y', 'gauss1');
        hwa(i) = hwhm(ff.c1);
        fits(i).axial.fit = ff;
        fits(i).axial.gof = gof;
        fits(i).axial.yvals = y;
        fits(i).axial.xvals = z;
        fits(i).axial.xvals_adj = z-ff.b1;
        fits(i).axial.hwhm = hwr(i);
    catch
        hwa(i) = nan;
        fits(i).axial = nan;
    end
end


%% do exclusions

% remove cells that didn't fit
excl = isnan(hwr) | isnan(hwa);

% remove cells that don't respond well
failures_r = max(rad(:,1:3),[],2) < stim_success;
failures_a = max(ax(:,1:3),[],2) < stim_success;
excl = excl | failures_r' | failures_a';

hwr = hwr(~excl);
hwa = hwa(~excl);
fits = fits(~excl);


%% plot them all initially

figure(21)
clf
hold on

jitter = 0.05;

% radial
yvals = hwr;
xvals = jitter*randn(numel(yvals),1) + 1;
s = scatter(xvals, yvals, 'filled');
s.MarkerEdgeColor = 'none';
s.MarkerFaceAlpha = 0.75;

% axial
yvals = hwa;
xvals = jitter*randn(numel(yvals),1) + 2;
s = scatter(xvals, yvals, 'filled');
s.MarkerEdgeColor = 'none';
s.MarkerFaceAlpha = 0.75;

xlim([0.25 2.75])

title('PPSF')
xticks([1 2])
xticklabels({'Radial', 'Axial'})
ylabel('HWHM')

%% view the actual cirve fits

one_by_one = 1;
center_fit = 1;


figure(22)
clf

for i=1:numel(hwa)
    % radial
    ff = fits(i).radial.fit;
    y = fits(i).radial.yvals;
    if center_fit
        x = fits(i).radial.xvals_adj;
    else
        x = fits(i).radial.xvals;
    end

    subplot(1,2,1)
    hold on

    plot(x,y, 'ok')
    plot(min(x):max(x),feval(ff,min(x):max(x)), 'r')
    title(['Radial, idx ' num2str(i)])

    % axial
    ff = fits(i).axial.fit;
    y = fits(i).axial.yvals;
    if center_fit
        z = fits(i).axial.xvals_adj;
    else
        z = fits(i).axial.xvals;
    end

    subplot(1,2,2)
    hold on

    plot(z,y, 'ok')
    plot(min(z):max(z),feval(ff,min(z):max(z)), 'r')
    title(['Axial, idx ' num2str(i)])

    disp(' ')
    disp(['Radial PPSF (hwhm): ' num2str(hwr(i))])
    disp(['Axial PPSF (hwhm):  ' num2str(hwa(i))])

    if one_by_one
        disp(['Cell idx # ' num2str(i)])
        pause
        clf
    end
end

disp(' ')
disp(['Mean radial PPSF (HWHM): ' num2str(mean(hwr))])
disp(['Median radial PPSF (HWHM): ' num2str(median(hwr))])
disp(['Mean axial PPSF (HWHM): ' num2str(mean(hwa))])
disp(['Median axial PPSF (HWHM): ' num2str(median(hwa))])
disp(' ')

%% correlate them
figure(23)
clf
scatter(hwr, hwa, 'filled', 'MarkerFaceAlpha', 0.75)
xlabel('Radial PPSF')
ylabel('Axial PPSF')
ylim([0 50])
xlim([0 50])
[rho, pval] = corr(hwr', hwa');

%% overall fit

center_fit = 1;
scaled = 1;
show_err = 1;

% radial
if center_fit
    xs = arrayfun(@(x) x.radial.xvals_adj, fits, 'UniformOutput', 0);
else
    xs = arrayfun(@(x) x.radial.xvals, fits, 'UniformOutput', 0);
end
xs = cell2mat(xs)';

ys = arrayfun(@(x) x.radial.yvals, fits, 'UniformOutput', 0);
if scaled
    ys = cellfun(@(x) x/max(x), ys, 'UniformOutput', 0);
end
yvals_r = cell2mat(ys');
ys = cell2mat(ys)';


rf = fit(xs, ys, 'gauss1');


% axial
if center_fit
    zs = arrayfun(@(x) x.axial.xvals_adj, fits, 'UniformOutput', 0);
else
    zs = arrayfun(@(x) x.axial.xvals, fits, 'UniformOutput', 0);
end
zs = cell2mat(zs)';

ys = arrayfun(@(x) x.axial.yvals, fits, 'UniformOutput', 0);
if scaled
    ys = cellfun(@(x) x/max(x), ys, 'UniformOutput', 0);
end
yvals_a = cell2mat(ys');
ys = cell2mat(ys)';

af = fit(zs, ys, 'gauss1');


% figure
figure(24)
clf

disp(' ')
disp('Overall fits....')

% RADIAL
subplot(1,2,1)
hold on

ppsf_rad_all = hwhm(rf.c1);

xr = -100:0.1:100;
plot(xr, feval(rf, xr), 'LineWidth',2)

x_int = feval(rf, ppsf_rad_all);
p = line([ppsf_rad_all ppsf_rad_all], [0 x_int]);
p.LineStyle = '--';

p = line([0 ppsf_rad_all], [x_int x_int]);
p.LineStyle = '--';

xline(0, '--')
xlim([-60 60])
title('Overall Radial Fit')

if show_err
    xplot = -3:3:30;
    if center_fit
        xplot = xplot-rf.b1;
    end
    
    m = mean(yvals_r,1);
    err = std(yvals_r,[],1)/(sqrt(size(yvals_r,2)));
    e = errorbar(xplot, m, err);
    e.LineWidth = 1.5;
    e.Marker = 'o';
    e.MarkerSize = 4;
    e.MarkerFaceColor = 'k';
    e.MarkerEdgeColor = 'none';
    e.LineStyle = 'none';
    e.Color = 'k';
end


disp(['PPSF (HWHM) Radial: ' num2str(ppsf_rad_all)])



% AXIAL
subplot(1,2,2)
hold on
ppsf_ax_all = hwhm(af.c1);

zr = -100:0.1:100;
plot(zr, feval(af, zr), 'LineWidth',2)

x_int = feval(af, ppsf_ax_all);
p = line([ppsf_ax_all ppsf_ax_all], [0 x_int]);
p.LineStyle = '--';

p = line([0 ppsf_ax_all], [x_int x_int]);
p.LineStyle = '--';

xline(0, '--')
xlim([-60 60])
title('Overall Axial Fit')

if show_err
    m = mean(yvals_a,1);
    err = std(yvals_a,[],1)/(sqrt(size(yvals_a,2)));
    xplot = -6:6:60;
    if center_fit
        xplot = xplot-af.b1;
    end
    e = errorbar(xplot, m, err);
    e.LineWidth = 1.5;
    e.Marker = 'o';
    e.MarkerSize = 4;
    e.MarkerFaceColor = 'k';
    e.MarkerEdgeColor = 'none';
    e.LineStyle = 'none';
    e.Color = 'k';
end

disp(['PPSF (HWHM) Axial: ' num2str(ppsf_ax_all)])

%% find lowest 10th percentile of activation

percentile = 0.1;

xr = 0:0.1:100;
vals = feval(rf, xr);
y_min = find(vals<percentile*max(vals),1);
rad_min = xr(y_min);
disp(['Radial distance for ' num2str(percentile) ' activation: ' num2str(rad_min)])

vals = feval(af, xr);
y_min = find(vals<percentile*max(vals),1);
ax_min = xr(y_min);
disp(['Axial distance for ' num2str(percentile) ' activation: ' num2str(ax_min)])

%% example cell

% current example is 13

ex = 13;

figure(25)
clf
% radial
subplot(1,2,1)
hold on
y = rad2(ex,:)';
x1 = xs(ex,:)';
ff = fit(x1,y, 'gauss1');
hwhm_rad = (2*sqrt(2*log(2))*ff.c1/sqrt(2))/2;
p1 = plot(x1,y,'o');
p2 = plot(min(x1):max(x1),feval(ff,min(x1):max(x1)));
p2.Color= p1.Color;
p2.LineWidth = 1;
title(['Radial ' num2str(hwhm_rad) ' \mum'])
ylabel('\DeltaF/F')
xlabel('Distance \mum')
xline(0, 'k--')
xlim([-3 33])

% axial
subplot(1,2,2)
hold on
y = ax2(ex,:)';
z1 = zs(ex,:)';
ff = fit(z1,y, 'gauss1');
hwhm_ax = (2*sqrt(2*log(2))*ff.c1/sqrt(2))/2;
p1 = plot(z1, y,'o');
p2 = plot(min(z1):max(z1),feval(ff,min(z1):max(z1)));
p2.Color= p1.Color;
p2.LineWidth = 1;
title(['Axial ' num2str(hwhm_ax) ' \mum'])
ylabel('\DeltaF/F')
xlabel('Distance \mum')
xline(0, 'k--')
xlim([-7 65])























