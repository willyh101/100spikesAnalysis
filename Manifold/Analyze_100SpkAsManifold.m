%% Load Experiments/Setup
clear
close all

masterTic = tic;

addpath(genpath('100spikesAnalysis'))
%% loadLists

ai203100spkLoadList
% oriLoadList;
% % allLoadList;

% loadPath = 'path/to/outfiles/directory';
loadPath = 'T:\Outfiles\caiman';

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
%% restrict Cells to use
opts.minMeanThreshold = -0.25;
opts.maxMeanThreshold = inf;

opts.verbose =0;
[All, cellExcludeResults] = cellExcluder(All,opts); 
allResults = cat(1,cellExcludeResults{:});
disp(['In total ' num2str(sum(allResults)) ' Cells Excluded. ' num2str(mean(allResults)*100,2) '%']);
disp(['Overall ' num2str(sum(~allResults)) ' Cells Passed!'])

opts.minNumCellsInd= 0;
tooFewCellsInds = cellfun(@(x) sum(~x)<opts.minNumCellsInd,cellExcludeResults);
disp([ num2str(sum(tooFewCellsInds)) ' inds have < ' num2str(opts.minNumCellsInd) ' cells, and should be exccluded']);

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
ensHNumber = outVars.ensHNumber; %in order of uniqueStims


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
opts.FractionMissable = 0.75;%0.33; %what percent of targets are missable before you exclude the ens
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
opts.IndsVisThreshold = 0.1;0.05; %default 0.05

highVisPercentInd = ~ismember(ensIndNumber,find(visPercent<opts.IndsVisThreshold)); %remove low vis responsive experiments
lowRunInds = ismember(ensIndNumber,find(percentLowRunTrials>0.5));
lowCellCount = ismember(ensIndNumber,find(tooFewCellsInds));


%exclude certain expression types:
uniqueExpressionTypes = outVars.uniqueExpressionTypes;
excludedTypes ={'AAV CamK2' 'neo-IV Tre 2s' 'IUE CAG' 'SepW1 CAG 2s' };


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
    ...& numMatchedTargets >= 3 ...
    ...& ensembleOneSecond ... %cuts off a lot of the earlier
    ...& numCellsEachEns==10 ...
    ... & ensDate > 210420 ...
    ...& outVars.hzEachEns == 10 ...
    ...& outVars.hzEachEns >= 9 & outVars.hzEachEns <= 12 ...
    & ~lowCellCount ...
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
disp(['Fraction of Ens high Cell Count: ' num2str(mean(~lowCellCount))]);


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

%% 
ensToUse = outVars.ensemblesToUse & outVars.ensOSI >0.7; 


%% quick determine stimmable

for ind =1:numel(All)
    
    stimTimes = All(ind).out.stm.holoRequest.bigListOfFirstStimTimes(:,1); 
    FR = All(ind).out.info.FR;
    
    dataToUse = All(ind).out.stm.zdfData; 
    
    tcs = All(ind).out.stm.targetedCells;
    for i = 1:numel(tcs)
        t = tcs(i);
        if ~isnan(t)
            st = round(stimTimes(i)*FR);
            
%             alignedDat = 
        end
    end
    
end

%% Zscore together

useCaiman = 1;

for ind=1:numel(All)
    if useCaiman
        dat1 = All(ind).out.vis.caiman_matched;
        dat2 = All(ind).out.exp.caiman_matched;
    else
        dat1 = All(ind).out.vis.zdfData;
        dat2 = All(ind).out.exp.zdfData;
    end
    
    sz1 = size(dat1);
    sz2 = size(dat2);
    mnFrames = min(sz1(2),sz2(2));
    
       fullDat = cat(3,dat1(:,1:mnFrames,:),dat2(:,1:mnFrames,:));

    
    if useCaiman
        sz = size(fullDat);
        fullDatUnRoll = reshape(fullDat,[sz(1) sz(2)*sz(3)]);
        
        fullzdfDat = zscore(fullDatUnRoll,[],2);
        fullzdfDat = reshape(fullzdfDat,sz);
    else
        [fulldfDat, fullzdfDat] = computeDFFwithMovingBaseline(fullDat);
    end
    
    
    dat1 = fullzdfDat(:,:,1:sz1(3));
    dat2 = fullzdfDat(:,:,sz1(3)+1:end);
    
    
    All(ind).out.anal.dat1 = dat1;
    All(ind).out.anal.dat2 = dat2;
    All(ind).out.anal.fullzdfDat = fullzdfDat; 
    
    
end

%% Plot Cosine Similarity 

numEns = numel(outVars.ensemblesToUse);

ensCrit = outVars.ensemblesToUse & outVars.ensOSI>0.7;

recWindow = 6:18;


useBothVisAndHoloForPCA=1;

EnsCosSim=[];
for i = 1:numEns
    if ensCrit(i)
        
        ind = outVars.ensIndNumber(i);
           
        htgs = unique([All(ind).out.exp.holoTargets{:}]);
        htgs(isnan(htgs))=[];
        cellsToUse = htgs;
        

        dat1 = All(ind).out.anal.dat1;
        dat2 = All(ind).out.anal.dat2; 
        
        try
            lowRunTrials = mean(All(ind).out.mani.runVal(:,6:12),2)<2;
            trialsToUseExp = All(ind).out.mani.lowMotionTrials & lowRunTrials';
        catch
            lowRunTrials = mean(All(ind).out.exp.runVal(:,6:12),2)<2;
            trialsToUseExp = All(ind).out.exp.lowMotionTrials & lowRunTrials';
        end
        lowRunTrialsVis = mean(All(ind).out.vis.runVal(:,6:12),2)<2;
        trialsToUseVis = All(ind).out.vis.lowMotionTrials & lowRunTrialsVis';
        
        visID = All(ind).out.vis.visID;
        
        visVector=zeros([numel(cellsToUse) 9]);
        for v=1:9
            visVector(:,v) = mean(mean(dat1(cellsToUse,recWindow, trialsToUseVis & visID ==v),2),3); 
        end
        
        s = outVars.ensHNumber(i);
        stimID = All(ind).out.exp.stimID; 
        us = unique(stimID);
        sID = us(s);
        
        
        holoVector = mean(mean(dat2(cellsToUse,recWindow, trialsToUseExp & stimID ==sID),2),3); 

        for v=1:9
            cosSim(v) = cosine_similarity(visVector(:,v),holoVector);
        end
        
        figure(1);clf
        ax(1) = subplot(1,3,1);
        imagesc(visVector)
        title('Vis')
        axis equal
        box off
        colorbar
        ax(2) = subplot(1,3,2);
        imagesc(holoVector)
        title('Holo')
        axis equal
        box off
        colorbar
        same_color_scale(ax)
        
       subplot(1,3,3)
       plot(cosSim)
       xticks(1:9)
       xticklabels({'grey' 0:45:315})
       xlabel('Oris')
       ylabel('Cosine Similarity')
       
       
        prefOri = outVars.ensPO(i);
        if isnan(prefOri)
            prefOri = 0;
        end
        
        orisPos = 0:45:315;
        
        %permute so preferred is in first position
%         p = find(orisPos==prefOri);
        [mx p] = max(cosSim(2:end));
        
        EnsCosSim(i,:) = [cosSim(1) circshift(cosSim(2:end),-1*(p-1))];
        
 
%     pause
        
        
    end
end

    figure(2);clf
    
    subplot(2,1,1)
    imagesc(EnsCosSim(ensCrit,:))
    caxis([-0.2 0.75])
    subplot(2,1,2)
    
    plot(mean(EnsCosSim(ensCrit,:)))
    
    EnsCosSimPreferred = EnsCosSim(ensCrit,2);
    EnsCosSimOther = EnsCosSim(ensCrit,[1 3:end]);

    disp(['Pref Ori Mean: ' num2str(mean(EnsCosSimPreferred)) ' +/- ' num2str(ste(EnsCosSimPreferred)) ])
    disp(['Other Ori Mean: ' num2str(mean(EnsCosSimOther(:))) ' +/- ' num2str(ste(EnsCosSimOther(:))) ])
    pPrefVOther = ranksum(EnsCosSimOther(:),EnsCosSimPreferred);
disp(['P Prev vs Other: ' num2str(pPrefVOther)])

    

    
    figure(4);clf
    plotSpread({EnsCosSimPreferred(:) EnsCosSimOther(:)})
    set(gca,'TickDir','out')
    %%
    
ensCrit = outVars.ensemblesToUse & outVars.ensOSI>0.7;% & outVars.numCellsEachEns == 3;
    
    useBothVisAndHoloForPCA=0;
    rangeToAnalyze =1:24; %for PCA
    plotRange = 1:12; %1:15 orig
    respRange = 6:18; %6:18 orig
    smoothFactor =3;  
    plotMerge = 1; 
    plotNoStim =0;
    
    plotLastHalfOnly =0;
    plotFirstHalfOnly=0; 
    
    plotExamples=0;
    
    pauseEachInd=1;
    
    for ind = IndsUsed
        %%
        visID = All(ind).out.vis.visID;
        uv = unique(All(ind).out.vis.visID);
        
        dat1 = All(ind).out.anal.dat1;
        dat2 = All(ind).out.anal.dat2;
        
        htgs = unique([All(ind).out.exp.holoTargets{:}]);
        htgs(isnan(htgs))=[];
        cellsToUse = htgs;
        numel(cellsToUse)
        
        dat1 = dat1(:,rangeToAnalyze,:);
        dat2 = dat2(:,rangeToAnalyze,:);
        %     dat1 = All(ind).out.vis.zdfData(:,1:18,:);
        %     dat2 = All(ind).out.exp.zdfData(:,1:18,:);
        
        sz1 = size(dat1);
        sz2 = size(dat2);
        mnFrames = min(sz1(2),sz2(2));
        
        try
            lowRunTrials = mean(All(ind).out.mani.runVal(:,6:12),2)<2;
            trialsToUseExp = All(ind).out.mani.lowMotionTrials & lowRunTrials';
        catch
            lowRunTrials = mean(All(ind).out.exp.runVal(:,6:12),2)<2;
            trialsToUseExp = All(ind).out.exp.lowMotionTrials & lowRunTrials';
        end
        lowRunTrialsVis = mean(All(ind).out.vis.runVal(:,6:12),2)<2;
        trialsToUseVis = All(ind).out.vis.lowMotionTrials & lowRunTrialsVis';
        
        dat1(:,:,~trialsToUseVis)=nan;
        dat2(:,:,~trialsToUseExp)=nan;
        %     datToUse = cat(3,dat1(:,1:mnFrames,:),dat2(:,1:mnFrames,:)); %All(ind).out.vis.allData;
        if useBothVisAndHoloForPCA
            datToUse = cat(3,dat1,dat2);%All(ind).out.vis.zdfData; %All(ind).out.vis.allData;
        else
            datToUse = dat1;
        end
        
        datToUse = datToUse(cellsToUse,:,:); %All(ind).out.vis.allData;
        bData = mean(datToUse(:,1:5,:),2);
        datToUse = datToUse-bData; %baseline
        %     datToUse = smoothdata(datToUse,2); %smooth
        
        datToUse=datToUse(:,6:end,:);
        
        sz = size(datToUse);
        datToUseMakePCA = reshape(datToUse,[sz(1) sz(2)*sz(3)]);
        datToUseMakePCA = smoothdata(datToUseMakePCA,2); %smooth
        
        [coeff,score,latent,tsquared,explained,mu] = pca(datToUseMakePCA');
        cumsum(explained);
        
        %data to make vis Epoch
        datToUse = dat1; %All(ind).out.vis.zdfData; %All(ind).out.vis.allData;
        
        datToUse = datToUse(cellsToUse,:,:); %All(ind).out.vis.allData;
        bData = mean(datToUse(:,1:5,:),2);
        datToUse = datToUse-bData; %baseline
        %     datToUse = smoothdata(datToUse,2); %smooth
        datToUse = movmean(datToUse,smoothFactor,2); %smooth
        
        sz = size(datToUse);
        datToUse2 = reshape(datToUse,[sz(1) sz(2)*sz(3)]);
        
        % score of vis Data
        visScore = coeff' * datToUse2;
        visScore = reshape(visScore,[sz(1) sz(2) sz(3)]);
        
        colors = colorMapPicker(5,'plasma');
        colors = colors(1:4);
        colors = cat(2,{rgb('grey')},colors,colors);
        
        
        %data to make holo epoch
        hdatToUse = dat2; % All(ind).out.exp.zdfData;
        
        hdatToUse = hdatToUse(cellsToUse,:,:); %All(ind).out.vis.allData;
        bhData = mean(hdatToUse(:,1:5,:),2);
        hdatToUse = hdatToUse-bhData; %baseline
        %     hdatToUse = smoothdata(hdatToUse,2); %smooth
        hdatToUse = movmean(hdatToUse,smoothFactor,2); %smooth
        
        
        sz = size(hdatToUse);
        hdatToUse2 = reshape(hdatToUse,[sz(1) sz(2)*sz(3)]);
        
        newHScore = coeff' * hdatToUse2;
        % newScore = reshape(score,[sz(2) sz(3) sz(1)]);
        newHScore = reshape(newHScore,[sz(1) sz(2) sz(3)]);
        
        %define Axis to plot
        ax1 = 1;
        ax2 = 2;
        ax3 = 3;
        
        
        figure(3);clf
        hold on
        if plotNoStim
            datToPlot = squeeze(mean(visScore(:,plotRange,visID==uv(1) & trialsToUseVis),3));
            p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
            p.Color = colors{1};
            p.LineWidth = 2;
        end
        
        %plot vis
        for i=1:9;%[1 3 5 7 9];%[1 8 6 4 2];%1:9;%numel(uv)
            
            datToPlot = squeeze(mean(visScore(:,plotRange,visID==uv(i)& trialsToUseVis),3));
            
            if ~plotMerge
                if (~plotLastHalfOnly && ~plotFirstHalfOnly) || (plotLastHalfOnly && i>=6) || (plotFirstHalfOnly && i<6)
                    sz = size(datToPlot);
                    %    datToPlot = smooth(datToPlot,3);
                    datToPlot = reshape(datToPlot,sz);
                    %         p = plot3(datToPlot(:,1),datToPlot(:,2),datToPlot(:,3));
                    p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                    p.Color = colors{i};
                    p.LineWidth = 2;
                end
                
            else
                if i>1 && i<6
                    datToPlot = squeeze(mean(visScore(:,plotRange,(visID==uv(i) | visID==uv(i+4)) & trialsToUseVis),3));
                    p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                    p.Color = colors{i};
                    p.LineWidth = 2;
                end
            end
            
            
        end
         
        
        stimID = All(ind).out.exp.stimID;       
        us = unique(stimID);
        
        
        for i=1:numel(outVars.ensemblesToUse)
            if outVars.ensemblesToUse(i) && outVars.ensIndNumber(i)==ind && ensCrit(i)
                s = outVars.ensHNumber(i); %i think this is the order of resp
                %                h = All(ind).out.exp.stimParams.roi{s};
%                 us = unique(All(ind).out.exp.stimID);
%                 sID = us(s);
                PO = outVars.ensPO(i); %   outVars.prefOris{ind}(s);
                if isnan(PO)
                    PO=1;
                else
                    PO = find(PO==0:45:315)+1;
                end
                datToPlot = squeeze(mean(newHScore(:,plotRange,stimID==us(s) & trialsToUseExp),3));
                p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                p.Color = colors{PO}; %rgb('red'); %colors{i};
                p.LineStyle = ':';
                p.LineWidth = 3;
                
            end
        end
            title(['Ind: ' num2str(ind)]);
    xlabel(['PC' num2str(ax1)])
    ylabel(['PC' num2str(ax2)])
    zlabel(['PC' num2str(ax3)])
    
    if pauseEachInd
        pause
    end
    end
    
    
    