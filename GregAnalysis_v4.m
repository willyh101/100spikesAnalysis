%% Start Here
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('GregAnalysis_v4.m')); 
addpath(genpath(folder));
rmpath(folder)

%%

% Identifies a subset of outfiles to load
oriLoadList;
loadPath = '/Users/gregoryhandy/Research_Local/outputdata1';

%% Load the outfiles ('Raw data')
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

%% Fix Known Errors
%CAUTION! ERRORS WILL OCCUR IF YOU RUN MORE THAN ONCE!
[All] = allLoadListErrorFixer(All,loadList);

%% clean Data, and create fields.

opts.FRDefault=6;
opts.recWinRange = [0.5 1.5];% %from vis Start in s [1.25 2.5];

%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull

%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded

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

%% Key function doing distance metric

% Added in some center of mass distance: e.g., All(30).out.anal.centerDistbyHolo
[All, outVars] = meanMatrixVisandCorr_GH(All,opts,outVars); %one of the main analysis functions

%%
visPercent = outVars.visPercent;
outVars.visPercentFromExp = visPercent;
ensIndNumber =outVars.ensIndNumber;


%% Optional: Calc pVisR from Visual Epoch [CAUTION: OVERWRITES PREVIOUS pVisR]
[All, outVars] = CalcPVisRFromVis(All,opts,outVars);
visPercent = outVars.visPercent;
outVars.visPercentFromVis = visPercent;

%% if there is a red section (will run even if not...)
[outVars] = detectShotRedCells(All,outVars);
ensHasRed = outVars.ensHasRed;

try
    arrayfun(@(x) sum(~isnan(x.out.red.RedCells)),All)
    arrayfun(@(x) mean(~isnan(x.out.red.RedCells)),All)
catch
end
%% Identify the Experiment type for comparison or exclusion
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;
ensExpressionType = indExpressionType(outVars.ensIndNumber);
outVars.ensExpressionType = ensExpressionType;

%% Missed Target Exclusion Criteria
%detects if too many targets were not detected in S2p
opts.FractionMissable = 0.33; %what percent of targets are missable before you exclude the ens
[outVars] = missedTargetDetector(All,outVars,opts);

ensMissedTargetF = outVars.ensMissedTargetF; %Fraction of targets per ensemble Missed
ensMissedTarget = outVars.ensMissedTarget; %Ensemble is unuseable
numMatchedTargets = outVars.numMatchedTargets;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ensemble Exclusion Section %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
opts.IndsVisThreshold = 0.20; %default 0.05

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<opts.IndsVisThreshold)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));

%exclude certain expression types:
uniqueExpressionTypes = outVars.uniqueExpressionTypes;
% excludedTypes = {'AAV CamK2'};
excludedTypes = {'AAV CamK2' 'Ai203'};
% excludedTypes = {'Ai203'};
% excludedTypes = {'AAV Tre'};
% excludedTypes = {};

exprTypeExclNum = find(ismember(uniqueExpressionTypes,excludedTypes));
excludeExpressionType = ismember(ensExpressionType,exprTypeExclNum);

% only include times where rate == numpulses aka the stim period is 1s.
ensembleOneSecond = outVars.numSpikesEachEns./outVars.numCellsEachEns == outVars.hzEachEns;

%spot to add additional Exclusions
excludeInds = ismember(ensIndNumber,[]); %Its possible that the visStimIDs got messed up

%Options
opts.numSpikeToUseRange = [75 150];[1 inf];[80 120];%[0 1001];
opts.ensStimScoreThreshold = 0.5; % default 0.5
opts.numTrialsPerEnsThreshold = 5; % changed from 10 by wh 4/23 for testing stuff

lowBaseLineTrialCount = ismember(ensIndNumber,find(numTrialsNoStimEns<opts.numTrialsPerEnsThreshold));

%% IMPORTANT: specifies the ensemble size to use
ensemblesToUse = numSpikesEachEns > opts.numSpikeToUseRange(1) ...
    & numSpikesEachEns < opts.numSpikeToUseRange(2) ...
    ...& highVisPercentInd ...
    & lowRunInds ...
    & ensStimScore > opts.ensStimScoreThreshold ... %so like we're excluding low success trials but if a holostim is chronically missed we shouldn't even use it
    & ~excludeInds ...
    & numTrialsPerEns > opts.numTrialsPerEnsThreshold ... ;%10;%&...
    & ~lowBaseLineTrialCount ...
    & ~ensHasRed ...
    & ~excludeExpressionType ...
    & ~ensMissedTarget ...
    & numMatchedTargets >= 3 ...
    ...& ensembleOneSecond ... %cuts off a lot of the earlier
    & numCellsEachEns==10 ...    
    ;
%& numCellsEachEns>10 ;

%%

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

disp(['Total Fraction of Ens Used: ' num2str(mean(ensemblesToUse))]);
disp([num2str(sum(ensemblesToUse)) ' Ensembles Included'])


%% Set Default Trials to Use
for ind=1:numExps
    trialsToUse = All(ind).out.exp.lowMotionTrials &...
        All(ind).out.exp.lowRunTrials &...
        All(ind).out.exp.stimSuccessTrial;% & ...
%         All(ind).out.exp.visID == 1 | ...
%         All(ind).out.exp.visID == 0;
    All(ind).out.anal.defaultTrialsToUse = trialsToUse;
end

%% Combine Ensemble Sizes 
numCellsEachEns = outVars.numCellsEachEnsBackup;
outVars.numCellsEachEns= numCellsEachEns;

SizeOption = 1; %0 true; 1 small medium and large; 2 all the same value

switch SizeOption
    case 1
        % Optional group ensembles into small medium and large
        numCellsEachEns = outVars.numCellsEachEnsBackup;
        numCellsEachEns(numCellsEachEns < 10) = 5;
        numCellsEachEns(numCellsEachEns > 10) = 20;
        
        outVars.numCellsEachEns= numCellsEachEns;
    case 2
        % set them all the same
        numCellsEachEns(numCellsEachEns >0 ) = 1;
        outVars.numCellsEachEns= numCellsEachEns;
end

%% Compare Distance responses
figure(102);clf

opts.distBins = 0:25:1000;

distTypes = {'min' 'geo' 'mean' 'harm' 'center of mass'};
for i =1:5
    disp(['working on ' distTypes{i}])
    opts.distType = distTypes{i}; %options: min geo mean harm
    CellToUseVar = [];
    [popRespDist] = popDistMaker_GH(opts,All,CellToUseVar,0);
    ax = subplot(2,3,i);
    opts.distAxisRange = [0 550]; %[0 350] is stand
    plotDistRespGeneric(popRespDist,outVars,opts,ax);
    title(distTypes{i})
end
disp('done')

%% Distance of Ensemble
[All, outVars] = defineDistanceTypes_GH(All, outVars);

%%
outVars.defaultColorMap = 'viridis';
plotEnsembleDistanceResponse_GH(outVars,100,1)

%%
clearvars dist_measures
dist_measures(1,:) = outVars.ensMeaD(outVars.ensemblesToUse);
% dist_measures(2,:) = outVars.ensGeoD(outVars.ensemblesToUse);
dist_measures(2,:) = outVars.ensMaxD(outVars.ensemblesToUse);
% dist_measures(4,:) = outVars.ensMinD(outVars.ensemblesToUse);
dist_measures(3,:) = outVars.ensCenterD(outVars.ensemblesToUse);
dist_measures(4,:) = outVars.ensCenterV(outVars.ensemblesToUse);

figure()
names = {'Mean','Max','Center D','Center Var'};
corrMatrix = corr(dist_measures').^2;
imagesc(corrMatrix)
colormap('gray')
colorbar
caxis([0 1])
set(gca,'fontsize',16)
xticks([1 2 3 4])
yticks([1 2 3 4])
set(gca,'XTickLabel',names)
set(gca,'YTickLabel',names)
title('R^2') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start of the OSI/Orientation Section %%%
%% Orientation Tuning and OSI %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI_GH(All, outVars);
[All, outVars] = calcTuningCircular_GH(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)
[All, outVars] = getEnsembleOSI_GH(All, outVars); % for ensembles specifically
%% plot vis things
opts.ensOSImethod = 'ensGOSI';% 'ensOSI'; 'meanEnsOSI'

plotOSIdists(outVars, opts);

%%
opts.distType = 'min';
[outVars] = grandDistanceMaker(opts,All,outVars);
numEns = numel(outVars.ensStimScore);

%%

clearvars osi_measures;
osi_measures(1,:) = outVars.ensOSI(outVars.ensemblesToUse);
osi_measures(2,:) = outVars.meanEnsOSI(outVars.ensemblesToUse);
osi_measures(3,:) = outVars.ensGOSI(outVars.ensemblesToUse);
osi_measures(4,:) = outVars.meanEnsGOSI(outVars.ensemblesToUse);
osi_measures(5,:) = 1-outVars.ensCircVarPO(outVars.ensemblesToUse);
bad_columns = find(sum(isnan(osi_measures))>0);
osi_measures(:,bad_columns) = [];

names = {'Ens OSI','Ave Ens OSI','Ens GOSI','Ave Ens GOSI','Ens 1-Circ Var PO'};

figure();
fancyPlotSpread(osi_measures',names);
set(gca,'fontsize',16)

figure()
corrMatrix = corr(osi_measures').^2;
imagesc(corrMatrix)
colormap('gray')
colorbar
caxis([0 1])
xticks([1 2 3 4 5])
yticks([1 2 3 4 5])
set(gca,'fontsize',16)
set(gca,'XTickLabel',names)
set(gca,'YTickLabel',names)
title('R^2') 


%% Plot Distance Plots by criteria
% More Criteria in this case set up as ensOSI, but its supposed to be
% flexible in case you want to change it later.

opts.distType = 'center of mass';
opts.distBins =[10:20:150];% [0:25:400];

%things to hold constant
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;%  &  outVars.ensOSI>0.75;;
criteria = outVars.ensOSI;outVars.ensMaxD;outVars.ensOSI; %outVars.ensMaxD;
useableCriteria = criteria(ensemblesToUse);
bins = [0 prctile(useableCriteria,25) prctile(useableCriteria,50) prctile(useableCriteria,75) max(useableCriteria)];

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

    cellToUseVar = ~outVars.offTargetRiskEns{i}...
        & outVars.pVisR{ind} < 0.1 ...
        & outVars.osi{ind} > 0.25 ...
        ;

        popToPlot(i,:) = popDistMakerSingle_GH(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
end
disp('Done')

figure(10);clf
figure(11);clf

diffsPossible = [0 45 90 135 180];
colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
for k = 1:numel(bins)-1
    figure(10);ax(k) =subplot(1,numel(bins)-1,k);
    title(['Ens Criteria: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ])
    popToPlotTemp = popToPlot;
    popToPlotTemp(criteria<bins(k) | criteria>bins(k+1),:)=NaN;
    [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(k));
    if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
        eHandle{1}.Color = colorListOri{k};
    end

    figure(11); tempax = subplot(1,1,1);
    [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
    eHandle{1}.Color = colorListOri{k};
end
linkaxes(ax)
figure(11);
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])
figure(10);
xlim([0 opts.distBins(end)])
ylim([-0.1 0.3])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Two Criteria Distance Plots
%%% Plot Distance Plots by Ensemble OSI and Distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.distType = 'min';
% opts.distType = 'center of mass';
opts.distBins =[15:20:150];% [0:25:400];

%adjust the Ensembles used 
ensemblesToUse = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10 ;%& outVars.meanEnsOSI>0.25;
% criteria = outVars.ensMaxD; criteria1_name{1} = 'Ens Max Dist'; bins = [0 400 500 inf];
% criteria = outVars.ensMeaD; criteria1_name{1} = 'Ens Mean Dist'; bins = [0 205.4 235.2 inf];
criteria = outVars.ensCenterD; criteria1_name{1} = 'Ens Center Dist'; bins = [0 140 180 inf];
% criteria = outVars.ensCenterV; criteria1_name{1} = 'Ens Center Var'; bins = [0 160 200 inf];
useableCriteria1 = criteria(ensemblesToUse);

% criteria2 = outVars.ensOSI; criteria2_name{1} = 'Ens OSI';
criteria2 = 1-outVars.ensCircVarPO; criteria2_name{1} = 'Ens PO 1-Circ Var';
useableCriteria2 = criteria2(ensemblesToUse);
bins2 = [0 0.5 inf]; 


figure(77)
subplot(2,1,1)
hold off
hist_data = histogram(useableCriteria1,20);
hold on
plot(bins(2)+[0:0.1:max(hist_data.Values)*1.1]*0,[0:0.1:max(hist_data.Values)*1.1],'k-','linewidth',2)
plot(bins(3)+[0:0.1:max(hist_data.Values)*1.1]*0,[0:0.1:max(hist_data.Values)*1.1],'k-','linewidth',2)

xlabel(criteria1_name{1})
set(gca,'fontsize',16)

subplot(2,1,2)
hold off
hist_data = histogram(useableCriteria2,20);
hold on
plot(bins2(2)+[0:0.1:max(hist_data.Values)*1.1]*0,[0:0.1:max(hist_data.Values)*1.1],'k-','linewidth',2)
plot(bins2(3)+[0:0.1:max(hist_data.Values)*1.1]*0,[0:0.1:max(hist_data.Values)*1.1],'k-','linewidth',2)

xlabel(criteria2_name{1})
set(gca,'fontsize',16)


%% Not entirely sure what this code does, but important for later code

ensIndNumber = outVars.ensIndNumber;
ensHNumber = outVars.ensHNumber;

clear popToPlot
for i=1:numEns %i know its slow, but All is big so don't parfor it
    if mod(i,round(numEns/10))==1
        fprintf('.')
    end

    if ensemblesToUse(i)
        ind = ensIndNumber(i);

        %adjust the responder cells included here
    cellToUseVar = ~outVars.offTargetRiskEns{i}...
         & outVars.pVisR{ind} < 0.05 ...
        ...& outVars.osi{ind} > 0.25 ...
        ... & outVars.posCellbyInd{i} ...
        ;

        popToPlot(i,:) = popDistMakerSingle(opts,All(ensIndNumber(i)),cellToUseVar,0,ensHNumber(i));

    else
        popToPlot(i,:) = nan([numel(opts.distBins)-1 1]);
    end
    
end
disp('Done')

%%
figure(13);clf

colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
c=0;
% closeData = -100+zeros(50,3);
total_count = 1;
for i = 1:numel(bins2)-1
    for k = 1:numel(bins)-1
        c=c+1;
        
        ax(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
        popToPlotTemp = popToPlot;
        if i == 1
        ensembleExcluder  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;
        else
        ensembleExcluder  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;    
        end
        popToPlotTemp(~ensembleExcluder,:)=NaN;
        
        [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(c));
             
        % New code
        nearToEnsData{total_count} = popToPlotTemp(~isnan(popToPlotTemp(:,1)),1);
        total_count = total_count + 1;
          
        figure(13)
        if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
            eHandle{1}.Color = colorListOri{k};
        end
        title({...
            %['Ens Dist: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            [criteria1_name{1} ': ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            %['Ens OSI: ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            [criteria2_name{1} ': ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]
            ['Num Ens: ' num2str(sum(ensembleExcluder & ensemblesToUse)) ] ...
            } )
        %ylabel('Pop Response (Mean \DeltaF/F)')
        ylabel('Pop Response')
        
        set(gca,'fontsize',16)
    end
    
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.2 0.25])


analyzeNearData_GH(nearToEnsData,bins2)


%%
figure(15);clf

colorListOri = colorMapPicker(numel(diffsPossible),'plasma');
clear ax
c=0;
for i = 1:1%numel(bins2)-1
    for k = 1:1%numel(bins)-1
        c=c+1;
        
%         ax(c) =subplot(numel(bins2)-1,numel(bins)-1,c);
%         ax(c) =subplot(1,numel(bins)-1,k);
        ax(c) =subplot(1,1,k);
        popToPlotTemp = popToPlot;
        ensembleExcluder  = criteria>=bins(k) & criteria<bins(k+1) & criteria2>=bins2(i) & criteria2<bins2(i+1) & ensemblesToUse;
        popToPlotTemp(~ensembleExcluder,:)=NaN;
        [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,ax(c));
        xlim([0 opts.distBins(end)])
        ylim([-0.2 0.25])
%         
%         figure(99)
        
        
        ensemblesToUse = outVars.ensemblesToUse;
        numCellsEachEns = outVars.numCellsEachEns;
        
        ensembles_size = unique(numCellsEachEns(ensemblesToUse));
        
        meanDat = mean(popToPlotTemp(find(numCellsEachEns==ensembles_size),:),'omitnan');
        
        semDat = std(popToPlotTemp(find(numCellsEachEns==ensembles_size),:),'omitnan')./...
            sqrt(sum(~isnan(popToPlotTemp(find(numCellsEachEns==ensembles_size),:))));
%         plot(opts.distBins(2:end)-length(opts.distBins)/2,mean(popToPlotTemp,'omitnan'),'.-','markersize',16)
        hold on;
        errorbar(opts.distBins(2:end)-(opts.distBins(2)-opts.distBins(1))/2,meanDat,semDat,'linewidth',2);
        xlim([0 opts.distBins(end)])
        ylim([-0.2 0.25])
        
        %%
        if numel(unique(outVars.numCellsEachEns(ensemblesToUse)))==1
            eHandle{1}.Color = colorListOri{k};
        end
        title({...
            ['Ens Dist: ' num2str(bins(k),2) ' to ' num2str(bins(k+1),2) ]...
            ['Ens OSI: ' num2str(bins2(i),2) ' to ' num2str(bins2(i+1),2) ]...
            ['Num Ens: ' num2str(sum(ensembleExcluder & ensemblesToUse)) ] ...
            } )
        
        ylabel('Pop Response (Mean \DeltaF/F)')
        
        %     figure(11); tempax = subplot(1,1,1);
        %     [eHandle] = plotDistRespGeneric(popToPlotTemp,outVars,opts,tempax);
        %     eHandle{1}.Color = colorListOri{k};
    end
end
linkaxes(ax)
xlim([0 opts.distBins(end)])
ylim([-0.2 0.25])

