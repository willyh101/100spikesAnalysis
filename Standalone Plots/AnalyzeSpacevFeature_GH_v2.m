%% Load Experiments/Setup
tic;

clear
close all

restoredefaultpath;
folder = fileparts(which('GregAnalysis_v4.m')); 
addpath(genpath(folder));
rmpath(folder)
% addpath(genpath('100spikesAnalysis'))
%% loadLists

oriLoadList;
% % allLoadList;

% loadPath = 'path/to/outfiles/directory';
% loadPath = 'T:\Outfiles';
loadPath = '/Users/gregoryhandy/Research_Local/outputdata1';

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

%% clean Data, and create fields.

opts.FRDefault=6;
opts.recWinRange = [0.5 1.5]; %[0.5 1.5];[1.5 2.5];%[0.5 1.5];% %from vis Start in s [1.25 2.5];


%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull
opts.stimSuccessByZ = 1; %if 0 calculates based on dataTouse, if 1 calculates based on zdfDat;


%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded
opts.runValPercent = 0.75; %percent of frames that need to be below run threshold

[All, outVars] = cleanData(All,opts);

ensStimScore        = outVars.ensStimScore;
hzEachEns           = outVars.hzEachEns;
numCellsEachEns     = outVars.numCellsEachEns;
numSpikesEachStim   = outVars.numSpikesEachStim;
percentLowRunTrials = outVars.percentLowRunTrials;
numSpikesEachEns    = outVars.numSpikesEachEns;
numSpikesEachCell   = outVars.numSpikesEachCell;

outVars.numCellsEachEnsBackup = outVars.numCellsEachEns;

names=[];
for Ind = 1:numel(All)
    names{Ind}=lower(strrep(All(Ind).out.info.mouse, '_', '.'));
end
outVars.names = names;

%% Make all dataPlots into matrixes of mean responses
%%Determine Vis Responsive and Process Correlation

opts.visAlpha = 0.05;

%oftarget risk params
opts.thisPlaneTolerance = 11.25;%7.5;%1FWHM%10; %in um;% pixels
opts.onePlaneTolerance = 22.5;%15;%2FWHM %20;
opts.distBins =  [0:25:1000]; [0:25:1000];
opts.skipVis =1;

[All, outVars] = meanMatrixVisandCorr(All,opts,outVars); %one of the main analysis functions

visPercent = outVars.visPercent;
outVars.visPercentFromExp = visPercent;
ensIndNumber =outVars.ensIndNumber;


%% REQUIRED: Calc pVisR from Visual Epoch [CAUTION: OVERWRITES PREVIOUS pVisR]
% Always do this!! not all experiments had full orientation data during the
% experiment epoch (but did during the vis epoch)
disp('Recalculating Vis stuff...')
opts.visRecWinRange = [0.5 1.5]; [0.5 1.5];
[All, outVars] = CalcPVisRFromVis(All,opts,outVars);
visPercent = outVars.visPercent;
outVars.visPercentFromVis = visPercent;

%% Misc Additional Variables:
%%RedSection: if there is a red section (will run even if not...)
[outVars] = detectShotRedCells(All,outVars);
ensHasRed = outVars.ensHasRed;

try
arrayfun(@(x) sum(~isnan(x.out.red.RedCells)),All)
arrayfun(@(x) mean(~isnan(x.out.red.RedCells)),All)
catch end
%%Identify the Experiment type for comparison or exclusion
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;
ensExpressionType = indExpressionType(outVars.ensIndNumber);
outVars.ensExpressionType = ensExpressionType;

%%Missed Target Exclusion Criteria
%detects if too many targets were not detected in S2p
opts.FractionMissable = 0.33; %what percent of targets are missable before you exclude the ens
[outVars] = missedTargetDetector(All,outVars,opts);

ensMissedTargetF = outVars.ensMissedTargetF; %Fraction of targets per ensemble Missed
ensMissedTarget = outVars.ensMissedTarget; %Ensemble is unuseable
numMatchedTargets = outVars.numMatchedTargets;
%%Determine Date
ensDate=[];
for i = 1:numel(ensIndNumber)
    ensDate(i) = str2num(All(ensIndNumber(i)).out.info.date);
end
outVars.ensDate=ensDate;

%%Identify duplicate holograms
[outVars] = identifyDuplicateHolos(All,outVars);

%% main Ensembles to Use section

numTrialsPerEns =[];numTrialsPerEnsTotal=[]; numTrialsNoStimEns=[];
for ind=1:numExps
    us=unique(All(ind).out.exp.stimID);

    for i=1:numel(us)
        trialsToUse = All(ind).out.exp.lowMotionTrials &...
            All(ind).out.exp.lowRunTrials &...
            All(ind).out.exp.stimSuccessTrial &...
            All(ind).out.exp.stimID == us(i) & ...
            (All(ind).out.exp.visID == 1 | All(ind).out.exp.visID == 0); %restrict just to no vis stim conditions

        numTrialsPerEns(end+1)=sum(trialsToUse);
        numTrialsPerEnsTotal(end+1) = sum(All(ind).out.exp.stimID == us(i));

        if i==1
            numTrialsNoStimEns(ind) = sum(trialsToUse);
        end
    end
end
numTrialsPerEns(numSpikesEachStim==0)=[];
numTrialsPerEnsTotal(numSpikesEachStim==0)=[];

%ID inds to be excluded
opts.IndsVisThreshold = 0.05; %default 0.05

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<opts.IndsVisThreshold)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));

%exclude certain expression types:
uniqueExpressionTypes = outVars.uniqueExpressionTypes;
excludedTypes ={'AAV CamK2' 'Ai203' 'neo-IV Tre 2s' 'IUE CAG' };


exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

% only include times where rate == numpulses aka the stim period is 1s.
ensembleOneSecond = outVars.numSpikesEachEns./outVars.numCellsEachEns == outVars.hzEachEns;

%spot to add additional Exclusions
excludeInds = ismember(ensIndNumber,[]); %Its possible that the visStimIDs got messed up
% excludeInds = ismember(ensIndNumber,[]); 

%Options
opts.numSpikeToUseRange = [90 110];[1 inf];[80 120];%[0 1001];
opts.ensStimScoreThreshold = 0.5; % default 0.5
opts.numTrialsPerEnsThreshold = 5; % changed from 10 by wh 4/23 for testing stuff

lowBaseLineTrialCount = ismember(ensIndNumber,find(numTrialsNoStimEns<opts.numTrialsPerEnsThreshold));


ensemblesToUse = numSpikesEachEns > opts.numSpikeToUseRange(1) ...
    & numSpikesEachEns < opts.numSpikeToUseRange(2) ...
    & highVisPercentInd ...
    & lowRunInds ...
    & ensStimScore > opts.ensStimScoreThreshold ... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    & ~excludeInds ...
    & numTrialsPerEns > opts.numTrialsPerEnsThreshold ... ;%10;%&...
    & ~lowBaseLineTrialCount ...
    & ~ensHasRed ...
    & ~excludeExpressionType ...
    & ~ensMissedTarget ...
    & numMatchedTargets >= 3 ...
    & ensembleOneSecond ... %cuts off a lot of the earlier
    & numCellsEachEns==10 ...
    ...& ensDate >= -210428 ...
    ...& outVars.hzEachEns == 10 ...
    ...& outVars.hzEachEns >= 9 & outVars.hzEachEns <= 12 ...
    ;

%%remove repeats
 [ensemblesToUse, outVars] = removeRepeatsFromEnsemblesToUse(ensemblesToUse,outVars);

indsSub = ensIndNumber(ensemblesToUse);
IndsUsed = unique(ensIndNumber(ensemblesToUse));

sum(ensemblesToUse)

outVars.ensemblesToUse      = ensemblesToUse;
outVars.IndsUsed            = IndsUsed;
outVars.indsSub             = indsSub;
outVars.numTrialsPerEns     = numTrialsPerEns;
outVars.highVisPercentInd    = highVisPercentInd;
outVars.lowRunInds           = lowRunInds;

%%Optional: Where are the losses comming from

disp(['Fraction of Ens correct Size: ' num2str(mean(numSpikesEachEns > opts.numSpikeToUseRange(1) & numSpikesEachEns < opts.numSpikeToUseRange(2)))]);
disp(['Fraction of Ens highVis: ' num2str(mean(highVisPercentInd))]);
disp(['Fraction of Ens lowRun: ' num2str(mean(lowRunInds))]);
disp(['Fraction of Ens high stimScore: ' num2str(mean(ensStimScore>opts.ensStimScoreThreshold))]);
disp(['Fraction of Ens high trial count: ' num2str(mean(numTrialsPerEns>opts.numTrialsPerEnsThreshold))]);
disp(['Fraction of Control Ens high trial count: ' num2str(mean(~lowBaseLineTrialCount))]);
disp(['Fraction of Ens No ''red'' cells shot: ' num2str(mean(~ensHasRed))]);
disp(['Fraction of Ens usable Expression Type: ' num2str(mean(~excludeExpressionType))]);
disp(['Fraction of Ens enough targets detected by s2p: ' num2str(mean(~ensMissedTarget))]);
disp(['Fraction of Ens number targets matched >=3: ' num2str(mean(numMatchedTargets >= 3))]);
disp(['Fraction of Ens Stim took 1s (aka correct stim Rate): ' num2str(mean(ensembleOneSecond))]);
disp(['Fraction of Ens that were not repeats: ' num2str(mean(~outVars.removedRepeats)) ]);

disp(['Total Fraction of Ens Used: ' num2str(mean(ensemblesToUse))]);
% disp([num2str(sum(ensemblesToUse)) ' Ensembles Included'])
disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);


%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials ...
        & All(ind).out.exp.lowRunTrials ...
        & All(ind).out.exp.stimSuccessTrial ...
        & (All(ind).out.exp.visID == 1 |  All(ind).out.exp.visID == 0 ) ...
            ;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end

%% Orientation Tuning and OSI

[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
[All, outVars] = getEnsembleOSI(All, outVars); % for ensembles specifically

%% Distance of Ensemble
[All, outVars] = defineDistanceTypes(All, outVars);
opts.distType = 'min';
[outVars] = grandDistanceMaker(opts,All,outVars);

%% Plot Space and Feature

opts.distType = 'min'; opts.distBins =[0:20:150];%[15:20:150];% [0:25:400];
% opts.distType = 'mean'; opts.distBins =[50:25:200];%[15:20:150];% [0:25:400];
opts.distAxisRange = [min(opts.distBins) max(opts.distBins)]; %[0 350];
opts.plotTraces =0;
plotSpaceAndFeature(All,outVars,opts)
%%
numEns = numel(outVars.ensStimScore);
ensemblesToUse = ...
    outVars.ensemblesToUse ...
    & outVars.numCellsEachEnsBackup==10 ...
    ...& outVars.meanEnsOSI>0.5 ...
    ;

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

    cellToUseVar = ~outVars.offTargetRiskEns{i}...
          & outVars.pVisR{ind} < 0.05 ...
          ... & outVars.osi{ind} > 0.25 ...
        ... & outVars.posCellbyInd{i} ...
        ;
    
        opts.distType = 'conn dist v2'; opts.distBins =[0:0.5:6 inf]; xlabel_name{1} = 'Ave # connections';
        [~,distValsConn{i},cellResp{i}] = popDistMakerSingle_GH_v3(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

        opts.distType = 'min'; opts.distBins =[25:10:500]; xlabel_name{1} = 'Min Dist to Ens';
        [~,distValsMin{i},~] = popDistMakerSingle_GH_v3(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
        
        opts.distType = 'mean'; opts.distBins =[100:50:500]; xlabel_name{1} = 'Mean Dist to Ens';
        [~,distValsMean{i},~] = popDistMakerSingle_GH_v3(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
        
        
        opts.distType = 'gauss'; opts.distBins =[100:50:500]; xlabel_name{1} = 'Gauss Dist to Ens';
        [~,distValsGauss{i},~] = popDistMakerSingle_GH_v3(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));
        
        
        ensSpread{i} = outVars.ensMaxD(i);
    else
        %popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
        distValsConn{i} = nan;
        distValsMin{i} = nan;
        distValsMean{i} = nan;
        cellResp{i} = nan;
        ensSpread{i} = nan;
    end
end
disp('Done')

%%

%%
conn_dist_all = [];
min_dist_all = [];
mean_dist_all = [];
gauss_dist_all = [];
cellResp_all = [];
ensSpread_all = [];
for ii = 1:length(ensemblesToUse)
    
    if ensemblesToUse(ii) == 1
        conn_dist_all = [conn_dist_all distValsConn{ii}];
        min_dist_all = [min_dist_all distValsMin{ii}];
        mean_dist_all = [mean_dist_all distValsMean{ii}];
        gauss_dist_all = [gauss_dist_all distValsGauss{ii}];
        cellResp_all = [cellResp_all cellResp{ii}'];
        ensSpread_all = [ensSpread_all ensSpread{ii}*ones(1,length(distValsConn{ii}))];
    end
    
end

%%

opts.distType = 'min'; opts.distBins =[15:10:300]; xlabel_name{1} = 'Min Dist to Ens';
figure(1234); clf; 
individual_cells_plot(min_dist_all,cellResp_all,opts,xlabel_name)

%%
%%
opts.distType = 'conn dist'; opts.distBins =[0:1:6]; xlabel_name{1} = 'Ave # Conn';
figure(123456); clf; 
individual_cells_plot(conn_dist_all,cellResp_all,opts,xlabel_name)

%%
opts.distType = 'gauss'; opts.distBins =[0:1:8]; xlabel_name{1} = 'Ave # Conn';
figure(123456); clf; 
individual_cells_plot(gauss_dist_all,cellResp_all,opts,xlabel_name)

%%
opts.distType = 'mean'; opts.distBins =[50:100:500]; xlabel_name{1} = 'Mean Dist to Ens';
figure(123456); clf; 
individual_cells_plot(mean_dist_all,cellResp_all,opts,xlabel_name)


%%
diffsPossible = [0 45 90 135 180];

criteria =  outVars.ensMaxD; outVars.ensMeaD;%   outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse); criteria1_name{1} ='Ens Max Dist';
% bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins = [0 400 500 inf];


% bins = [linspace(min(useableCriteria),max(useableCriteria),5)];
criteria2 =   outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria2(ensemblesToUse); criteria2_name{1} ='Ens OSI';
% bins2 = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];
bins2 = [0 0.3 0.7 inf];


% opts.distType = 'min'; opts.distBins =[0:25:300]; xlabel_name{1} = 'Min Dist to Ens';
% opts.distType = 'mean'; opts.distBins =[0:20:300]; xlabel_name{1} = 'Mean Dist to Ens';
opts.distType = 'conn dist v2'; opts.distBins =[0:1:6]; xlabel_name{1} = 'Ave # connections';
% opts.distType = 'gauss'; opts.distBins =[0:0.5:8]; xlabel_name{1} = 'Gaussian';
figure(133); clf;

colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
dataForStats =[];
c=0;
for i = 1:numel(bins2)-1
    for k = 1:numel(bins)-1
        figure(135);  hold on
        c=c+1;
        
        ax(c) =subplot(numel(bins2)-1,numel(bins)-1,c); hold on;
        ensembleSelecter  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;

        title({...
            [criteria1_name{1} ': ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            [criteria2_name{1} ': ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            ['Num Ens: ' num2str(sum(ensembleSelecter & ensemblesToUse)) ] ...
            } )
        ylabel('Pop Response (Mean \DeltaF/F)')
        set(gca,'fontsize',16)
        xlabel(xlabel_name)
        xlim([0 opts.distBins(end)])
        
        distValsMin_small{c} = [];
        cellResp_small{c} = [];
        for ii = 1:length(ensembleSelecter)
            if ensembleSelecter(ii) == 1
                if strcmp(opts.distType, 'mean') == 1
                    distValsMin_small{c} = [distValsMin_small{c} distValsMean{ii}];
                    cellResp_small{c} = [cellResp_small{c} cellResp{ii}'];
                    scatter(distValsMean{ii},cellResp{ii}','filled','k',...
                        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
                elseif strcmp(opts.distType, 'min') == 1
                    distValsMin_small{c} = [distValsMin_small{c} distValsMin{ii}];
                    cellResp_small{c} = [cellResp_small{c} cellResp{ii}'];
                    plot(distValsMin{ii},cellResp{ii}','.','color',[0 0.4470 0.7410]);
                elseif strcmp(opts.distType, 'gauss') == 1
                    distValsMin_small{c} = [distValsMin_small{c} distValsGauss{ii}];
                    cellResp_small{c} = [cellResp_small{c} cellResp{ii}'];
                    scatter(distValsGauss{ii},cellResp{ii}',[],distValsConn{ii},'filled',...
                        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
                else
                    distValsMin_small{c} = [distValsMin_small{c} distValsConn{ii}];
                    cellResp_small{c} = [cellResp_small{c} cellResp{ii}'];
                    scatter(distValsConn{ii},cellResp{ii}',[],'filled','k',...
                        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
                end  
            end
        end
        ylim_temp = ylim;
        plot([0 opts.distBins], [0 opts.distBins*0], 'k--', 'LineWidth',2)
        
        
        [y_avg, y_std, y_sem] = calAvg(distValsMin_small{c}, cellResp_small{c}, opts.distBins, diff(opts.distBins(1:2)));
    
        
        figure(133); hold on;
        ax2(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
        temp_nan = isnan(y_avg);
        if temp_nan(1)==1
            start_plot = 1+find(isnan(y_avg),1,'last');
            end_plot = length(y_avg);
        elseif temp_nan(end)==1
            start_plot = 1;
            end_plot = find(isnan(y_avg),1,'first')-1;
        else
            start_plot = 1;
            end_plot = length(y_avg);
        end
        plot(opts.distBins(start_plot:end_plot), y_avg(start_plot:end_plot), 'r-', 'LineWidth',2)
        plot([0 opts.distBins(start_plot:end_plot)], [0 opts.distBins(start_plot:end_plot)*0], 'k--', 'LineWidth',2)
        patch([opts.distBins(start_plot:end_plot) fliplr(opts.distBins(start_plot:end_plot))],...
            [y_avg(start_plot:end_plot)-y_sem(start_plot:end_plot) fliplr(y_avg(start_plot:end_plot)+y_sem(start_plot:end_plot))], 'r','FaceAlpha',.3)
        set(gca,'fontsize',16)
        xlabel(xlabel_name)
   
        xlim([0 opts.distBins(end)])
        
        
        figure(13433); hold on;
        ax3(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
        boxplot(cellResp_small{c}(distValsMin_small{c}<25))
        title(sprintf('Median %0.2f',median(cellResp_small{c}(distValsMin_small{c}<25))))
    end
end
figure(133);
linkaxes(ax2)
% xlim([0 opts.distBins(end)])
% xlim([0 50])
%%
figure(135);
linkaxes(ax)
xlim([0 50])

%%

[H, p] = ttest2((cellResp_small{1}(distValsMin_small{1}<25)),...
    cellResp_small{7}(distValsMin_small{7}<25));

[H, p] = ttest2((cellResp_small{9}(distValsMin_small{9}<25)),...
    cellResp_small{7}(distValsMin_small{7}<25));




%%

function [y_avg, y_std, y_sem] = calAvg(x, y, x_avg, dx, num_th)
% if nargin == 4; num_th = 20; end
if nargin == 4; num_th = 1; end

    y_avg = nan*zeros(size(x_avg));
    y_std = nan*zeros(size(x_avg));
    y_sem = nan*zeros(size(x_avg));
    
    for i = 1:length(x_avg)
        xids = logical((x>(x_avg(i)-dx/2)) .* (x<(x_avg(i)+dx/2)));
        if sum(xids(:)) >= num_th
            y_avg(i) = nanmean(y(xids));
            y_std(i) = nanstd(y(xids));
            y_sem(i) = nanstd(y(xids)) / sqrt(length(xids));
        end
    end
end
