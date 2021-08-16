%%for NOT ai203 mice
disp('   ')
pvis_cutoff = 0.2;

numExpts = numel(All);

clear osis visr percent_vis
c=0;

uniqueExpressionTypes = outVars.uniqueExpressionTypes;
% excludedTypes = {'AAV CamK2' 'Ai203'};
% excludedTypes = {'Ai203'};
% excludedTypes = {'Ai203'};
excludedTypes = {'AAV Tre'};
% excludedTypes = {};

exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

inds2use = unique(ensIndNumber(excludeExpressionType));
disp(['Total number of Ai203 FOVs ' num2str(numel(inds2use))])

for ind = inds2use
    out = All(ind).out;
    visPercent = out.anal.visPercent;
    
    if visPercent < pvis_cutoff
        continue
    end
    
    c = c + 1;
    osis{c} = out.anal.osi;
    visr{c} = out.anal.pVisR < 0.05;
    percent_vis(c) = visPercent;
end

osis = cell2mat(osis);
visr = cell2mat(visr);
osis = osis(~isnan(osis));
visr = visr(~isnan(osis));

osis = osis(visr);
not_ai203_osis = osis;
virus_vis_percent = percent_vis;

disp(['Total number of cells included: ' num2str(numel(osis))])

figure(244)
clf
subplot(2,2,1)
histogram(osis,20)
title(['Not Ai203, mean= ' num2str(mean(osis))])
ylabel('Number of cells')
xlabel('OSI')

disp(['Mean OSI: ' num2str(mean(osis))])
disp(['Median OSI: ' num2str(median(osis))])

%% for ai203 mice


numExpts = numel(All);

clear osis visr percent_vis
c=0;

uniqueExpressionTypes = outVars.uniqueExpressionTypes;
excludedTypes = {'Ai203'};
% excludedTypes = {'Ai203'};
% excludedTypes = {'Ai203'};
% excludedTypes = {'AAV Tre'};
% excludedTypes = {};

exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

inds2use = unique(ensIndNumber(excludeExpressionType));
disp(['Total number of virus FOVs ' num2str(numel(inds2use))])

for ind = inds2use
    out = All(ind).out;
    visPercent = out.anal.visPercent;
    
    if visPercent < pvis_cutoff
        continue
    end
    
    c = c + 1;
    osis{c} = out.anal.osi;
    visr{c} = out.anal.pVisR < 0.05;
    percent_vis(c) = visPercent;
    out.info
end

osis = cell2mat(osis);
visr = cell2mat(visr);
osis = osis(~isnan(osis));
visr = visr(~isnan(osis));

osis = osis(visr);
ai203_osis = osis;
ai203_vis_percent = percent_vis;

disp(['Total number of cells included: ' num2str(numel(osis))])

subplot(2,2,2)
histogram(osis,20)
title('Ai203')
ylabel('Number of cells')
xlabel('OSI')
title(['Ai203, mean= ' num2str(mean(osis))])

disp(['Mean OSI: ' num2str(mean(osis))])
disp(['Median OSI: ' num2str(median(osis))])

%%choose a random example expt and plot good vs bad OSI

curves2plot = horzcat(outVars.tuningCurves{:});
curves2plotSEM = horzcat(outVars.tuningCurvesSEM{:});
allOsis2use = horzcat(outVars.osi{:});
curves2plot = curves2plot(2:end,:);
curves2plotSEM = curves2plotSEM(2:end,:);
curves2plot = curves2plot - min(curves2plot);

bad2plot = curves2plot(:,allOsis2use < 0.3);
good2plot = curves2plot(:,allOsis2use > 0.6);

bad2plotSEM = curves2plotSEM(:,allOsis2use < 0.3);
good2plotSEM = curves2plotSEM(:,allOsis2use > 0.6);

badOSIs = allOsis2use(allOsis2use < 0.3);
goodOSIs = allOsis2use(allOsis2use > 0.6);

ex = randperm(size(bad2plot, 2), 1);
ex2plot = bad2plot(:, ex);
ex2plotSEM = bad2plotSEM(:, ex);
badOSIused = badOSIs(ex);

exgood = randperm(size(good2plot, 2), 1);
exgood2plot = good2plot(:, exgood);
exgood2plotSEM = good2plotSEM(:, exgood);
goodOSIused = goodOSIs(exgood);

oris = 0:45:315;

subplot(2,2,3)
e = errorbar(oris, ex2plot, ex2plotSEM);
e.Color = rgb('grey');
e.LineWidth = 1;
title(['Low OSI: ' num2str(badOSIused)])
xlabel('Orientation')
ylabel('z-scored dF/F')
xticks(oris)

subplot(2,2,4)
e = errorbar(oris, exgood2plot, exgood2plotSEM);
e.Color = rgb('grey');
e.LineWidth = 1;
title(['High OSI: ' num2str(goodOSIused)])
xlabel('Orientation')
ylabel('z-scored dF/F')
xticks(oris)

%%statistics for OSI
disp('.........')
[h,p] = ttest2(not_ai203_osis, ai203_osis);
% [p,h] = ranksum(not_ai203_osis, ai203_osis);
disp(['pVal all samples: ' num2str(p)])

notai203_sub = randsample(not_ai203_osis, 100);
yesai203_sub = randsample(ai203_osis,100);

[h,p] = ttest2(notai203_sub, yesai203_sub);
% [p,h] = ranksum(not_ai203_osis, ai203_osis);
disp(['pVal subsampled: ' num2str(p)])

mpval = zeros(1,10000);
for r=1:10000
    notai203_sub = randsample(not_ai203_osis, 100);
    yesai203_sub = randsample(ai203_osis,100);
    [h,p] = ttest2(notai203_sub, yesai203_sub);
    mpval(r) = p;
end
p = mean(mpval);
disp(['pVal jackknife 10,000x: ' num2str(p)])

ncells_ai203 = numel(ai203_osis);
not_ai203_ss = randsample(not_ai203_osis,ncells_ai203);
[h,p] = ttest2(ai203_osis, not_ai203_ss);
disp(['pVal subsample not ai203 only: ' num2str(p)])

mpval = zeros(1,10000);
for r=1:10000
    not_ai203_ss = randsample(not_ai203_osis,ncells_ai203);
    [h,p] = ttest2(not_ai203_ss, ai203_osis);
    mpval(r) = p;
end
disp(['pVal jackknife 10,000x sub not ai203 only: ' num2str(p)])

% ranksum aka mann-whitney u test
p = ranksum(not_ai203_osis, ai203_osis);
disp(['pVal mann-whitney u-test (aka ranksum): ' num2str(p)])

disp('........')

%%visual responses

figure(245)
clf

% for different FOVs
beeSwarmPlot({ai203_vis_percent, virus_vis_percent}, {'Ai203', 'AAV'});
[h,p] = ttest2(ai203_vis_percent, virus_vis_percent);
disp(['pVal ttest vis responsive percent: ' num2str(p)])







disp('.......')
disp(' ')