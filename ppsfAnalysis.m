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

%% error fixer
%CAUTION! ERRORS WILL OCCUR IF YOU RUN MORE THAN ONCE!
[All] = allLoadListErrorFixer(All,loadList);

%% Set Data To use
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

resp = nan([4,24,10]);


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
            temp = nan(1,10);
            temp(1:length(cellVal)) = cellVal;
            resp(c,i,:) = temp;
        end
    end
end

%% make into radial and axial trials
rshape = reshape(permute(resp, [1 3 2]), 40, 24);
rad = rshape(:,1:12);
ax = rshape(:,13:24);
%% fit and find hwhm

% resp is 4 x 24 x 10, holos, trials, cells
pxPerMu = 800/512;
% optoPerZ = 60/55;
x = [-3:3:30] * pxPerMu;
z = [-6:6:60];


figure(4)
clf


% plot x offsets
subplot(2,2,1)
hold on
title('Radial PPSF')
clear ppsf_rad err_rad off_r
k=0;
for i=1:size(rad,1)
    k=k+1;
    y = rad(i,:);
    if all(isnan(y))
        ppsf_rad(k) = nan;
        off_r(k) = nan;
        continue
    end
    try
        [ff, gof] = fit(z',y', 'gauss1');
        fwhm = 2*sqrt(2*log(2))*ff.c1/sqrt(2);
        hwhm = fwhm/2;   
        off_r(k) = ff.b1;
    catch
        hwhm = nan;
        off_r(k) = nan;
    end

    ppsf_rad(k) = hwhm;
    

    plot(x,y, 'ok')
    plot(min(x):max(x),feval(ff,min(x):max(x)), 'r')

    %         pause
end


% plot axial offsets
subplot(2,2,2)
hold on
title('Axial PPSF')
k=0;
clear ppsf_ax err_ax off_a
for i=1:size(ax,1)
    k=k+1;
    y = ax(i,:);
    if all(isnan(y))
        ppsf_ax(k) = nan;
        off_a(k) = nan;
        continue
    end
    
    try
        [ff, gof] = fit(z',y', 'gauss1');
        fwhm = 2*sqrt(2*log(2))*ff.c1/sqrt(2);
        hwhm = fwhm/2;    
        off_a(k) = ff.b1;
    catch
        hwhm = nan;
        ff.b1 = nan;
        off_a(k) = nan;
    end

    ppsf_ax(k) = hwhm;
    
    plot(z,y, 'ok')
    plot(min(z):max(z),feval(ff,min(z):max(z)), 'r')

end

disp(' ')
disp(['Mean radial PPSF (HWHM): ' num2str(nanmean(ppsf_rad))])
disp(['Median radial PPSF (HWHM): ' num2str(nanmedian(ppsf_rad))])
disp(['Mean axial PPSF (HWHM): ' num2str(nanmean(ppsf_ax))])
disp(['Median axial PPSF (HWHM): ' num2str(nanmedian(ppsf_ax))])
disp(' ')


subplot(2,2,3)
% r2_rad = arrayfun(@(x) x.rsquare, err_rad, 'UniformOutput', false);
% hist(ppsf_rad)
% hist(ppsf_rad(r2_rad>0.8))
plot(ppsf_rad, 'o')

subplot(2,2,4)
% r2_ax = arrayfun(@(x) x.rsquare, err_ax, 'UniformOutput', false);
% hist(ppsf_ax)
% hist(ppsf_ax(r2_ax>0.8))
plot(ppsf_ax, 'o')

 %% problem points
% pp = find(ppsf_rad > 50 | ppsf_ax > 50);
% pp = find(ppsf_rad < 50 | ppsf_ax < 50);
% 
% figure(5)
% clf
% 
% 
% for i=1:numel(pp)
%     p = pp(i);
%     subplot(1,2,1)
%     hold on
%     title('Radial PPSF')
%     y = rad(p,:);
%     p1 = plot(x, y, '-o');
%     [ff, gof] = fit(x',y', 'gauss1');
%     p2 = plot(min(x):max(x),feval(ff,min(x):max(x)));
%     p2.Color= p1.Color;
%     p2.LineWidth = 2;
% 
%     disp(['Radial.... ' num2str(ppsf_rad(p))])
%     gof
%     
%     subplot(1,2,2)
%     hold on
%     title('Axial PPSF')
%     y = ax(p,:);
%     p1 = plot(z, y, '-o');
%     [ff, gof] = fit(z',y', 'gauss1');
%     p2 = plot(min(z):max(z),feval(ff,min(z):max(z)));
%     p2.Color= p1.Color;
%     p2.LineWidth = 2;
% 
%     disp(['Axial.... ' num2str(ppsf_ax(p))])
%     gof
% 
% %     pause
% end

%% correlate ppsf

figure(6)
clf
scatter(ppsf_rad, ppsf_ax)
xlabel('Radial PPSF')
ylabel('Axial PPSF')

%% fit all and plot

do_zscore = 0;
stimSuccess = 0.25;

stimFailCells = rad(:,1) < stimSuccess;

excl = zeros(size(ppsf_ax));
excl = excl | stimFailCells';
excl = excl | ppsf_rad > 50 | ppsf_ax > 50;
% excl = ppsf_rad > 50 | ppsf_ax > 50;
% excl = ppsf_rad < 50 | ppsf_ax < 50;


xs = repmat(x,40,1);
xs = xs-off_r';
rad2 = rad(any(~isnan(xs),2) & ~excl',:);
xs = xs(any(~isnan(xs),2) & ~excl',:);
xs = reshape(xs',[],1);

if do_zscore
    rad2 = zscore(rad2);
end
rad2 = reshape(rad2',[],1);

rf = fit(xs, rad2, 'gauss1');




zs = repmat(z,40,1);
zs = zs-off_a';
ax2 = ax(any(~isnan(zs),2) & ~excl',:);
zs = zs(any(~isnan(zs),2) & ~excl',:);
zs = reshape(zs',[],1);

if do_zscore
    ax2 = zscore(ax2,[]);
end
ax2 = reshape(ax2',[],1);



af = fit(zs, ax2, 'gauss1');


figure(6)
clf

disp(' ')
disp('Overall fits....')

subplot(1,2,1)
hold on
ppsf_rad_all = (2*sqrt(2*log(2))*rf.c1/sqrt(2))/2;
title(['Overall Radial Fit ' num2str(ppsf_rad_all) '\mum'])
xr = min(x):max(x);
plot(xr, feval(rf, xr), 'LineWidth',2)
% scatter(xs,rad2)
disp(['PPSF (HWHM) Radial: ' num2str(ppsf_rad_all)])

subplot(1,2,2)
hold on
ppsf_ax_all = (2*sqrt(2*log(2))*af.c1/sqrt(2))/2;
title(['Overall Axial Fit ' num2str(ppsf_ax_all) '\mum'])
zr = min(z):max(z);
plot(zr, feval(af, zr), 'LineWidth',2)
% scatter(zs, ax2)
disp(['PPSF (HWHM) Axial: ' num2str(ppsf_ax_all)])

%% fit indiv and plot

do_zscore = 0;
stimSuccess = 0.25;

stimFailCells = rad(:,1) < stimSuccess;

excl = zeros(size(ppsf_ax));
excl = excl | stimFailCells';
excl = excl | ppsf_rad > 50 | ppsf_ax > 50;

% excl = zeros(size(ppsf_ax));
% excl = ppsf_rad > 50 | ppsf_ax > 50;
% excl = ppsf_rad > 500 | ppsf_ax > 500;


xs = repmat(x,40,1);
xs = xs-off_r';
rad2 = rad(any(~isnan(xs),2) & ~excl',:);
xs = xs(any(~isnan(xs),2) & ~excl',:);


figure(8)
clf

subplot(1,2,1)
hold on
clear hw_x
title('Radial')
for i=1:size(rad2,1)
    y = rad2(i,:)';
    if do_zscore
        y = zscore(y);
    end
    x1 = xs(i,:)';
    ff = fit(x1,y, 'gauss1');
    hw_x(i) = (2*sqrt(2*log(2))*ff.c1/sqrt(2))/2;
    p1 = plot(x1,y,'o');
    p2 = plot(min(x1):max(x1),feval(ff,min(x1):max(x1)));
    p2.Color= p1.Color;
    p2.LineWidth = 2;

end


zs = repmat(z,40,1);
zs = zs-off_a';
ax2 = ax(any(~isnan(zs),2) & ~excl',:);
zs = zs(any(~isnan(zs),2) & ~excl',:);


subplot(1,2,2)
hold on
clear hw_z
title('Axial')
for i=1:size(ax2,1)
    y = ax2(i,:)';
     if do_zscore
        y = zscore(y);
     end
    z1 = zs(i,:)';
    ff = fit(z1,y, 'gauss1');
    hw_z(i) = (2*sqrt(2*log(2))*ff.c1/sqrt(2))/2;
    p1 = plot(z1, y,'o');
    p2 = plot(min(z1):max(z1),feval(ff,min(z1):max(z1)));
    p2.Color= p1.Color;
    p2.LineWidth = 2;
end

disp(' ')
disp(['Radial ppsf (hwhm): ' num2str(mean(hw_x))])
disp(['Axial  ppsf (hwhm): ' num2str(mean(hw_z))])











