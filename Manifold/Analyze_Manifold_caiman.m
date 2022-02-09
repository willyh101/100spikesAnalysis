clear
close all

masterTic = tic;

addpath(genpath('100spikesAnalysis'))
%% loadLists

ai203ManifoldLoadListCaiman;
% manifoldLoadList;
% % allLoadList;

% loadPath = 'path/to/outfiles/directory';
% loadPath = 'T:\Outfiles';
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

%% create dfDatas if needed
numExps = numel(All);

for ind = 1:numExps
    if ~isfield(All(ind).out.vis,'dfData') | ~isfield(All(ind).out.vis,'zdfData')
        [dfData, zdfData]  = computeDFFwithMovingBaseline(All(ind).out.vis.allData);
        All(ind).out.vis.dfData = dfData;
        All(ind).out.vis.zdfData = zdfData;
    end
end

%% Set Data To Use

for ind = 1:numExps
%     
%   All(ind).out.exp.dataToUse = All(ind).out.exp.zdfData;
%     All(ind).out.vis.dataToUse = All(ind).out.vis.zdfData; 
%     All(ind).out.exp.dataToUse = All(ind).out.exp.caiman_matched;
%     All(ind).out.vis.dataToUse = All(ind).out.vis.caiman_matched;  
  All(ind).out.exp.dataToUse = All(ind).out.exp.caiman;
    All(ind).out.vis.dataToUse = All(ind).out.vis.caiman;  
end
%%


opts.FRDefault=6;
opts.recWinRange = [0.5 1.5]; %[0.5 1.5];[1.5 2.5];%[0.5 1.5];% %from vis Start in s [1.25 2.5];
opts.visRecWinRange = [0.5 1.5];

%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull
opts.stimSuccessByZ = 1; %if 0 calculates based on dataTouse, if 1 calculates based on zdfDat;


%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded
opts.runValPercent = 0.75; %percent of frames that need to be below run threshold

All = cleanDataVisOnly(All,opts);

%%
disp('Calculating Vis stuff...')
opts.visAlpha = 0.05;
[All, outVars] = CalcPVisRFromVis(All,opts);
visPercent = outVars.visPercent;
outVars.visPercentFromVis = visPercent;

disp('ready')

%% Power Calibration Sumamry
plotExamples = 0; 

tic;
satPowerForFitCellsTemp ={};
allPBal ={};
for ind = 1:numExps; %comment out the plot example part to run it as parfor
    disp(['Starting ind: ' num2str(ind)]);
    stm = All(ind).out.stm;
    FR1 = All(ind).out.info.FR;
    
    try
        stime = stm.holoRequest.bigListOfFirstStimTimes(:,1);
    catch
        stime = stm.holoRequest2.bigListOfFirstStimTimes(:,1);
    end
    stFr = round(stime*FR1);
    
    dataToUse1 = stm.caiman; %_matched;    stm.zdfData;
    sz = size(dataToUse1);
    
    trialsToUse = stm.lowMotionTrials & mean(stm.runVal,2)'<6;
    dataToUse1(:,:,~trialsToUse)=nan;
    
    tc =stm.targetedCells; 
    tc = (1:numel(stm.targetedCells))';
    
    calibCells  = ~isnan(tc) & ~isnan(stime);
    datAlign = nan([sum(calibCells) 12 sz(3)]);
    counter =1;
    for i=1:numel(tc);
        if calibCells(i);
            c = tc(i);
            st = stFr(i);
            if st+6 > sz(2)
                padSize = st+6-sz(2);
                
                if padSize > 12
                    datAlign(counter,:,:) = nan([1 12 sz(3)]);
                else
                    datpart = dataToUse1(c,st-5:end,:);
                    datAlign(counter,:,:) = padarray(datpart,[0 padSize 0],nan,'post');
                end
            else
                datAlign(counter,:,:) = dataToUse1(c,st-5:st+6,:);
            end
            counter=counter+1;
        end
    end
    
    powerList = stm.holoRequest.holoStimParams.powerList{:};
    outputsOrder = stm.outputsInfo.OutputOrder;
    stimID1 = stm.stimID;
    %     stimID = stimID(1:77);
    %     us = unique(stimID1);
    
    powers=[];
    for i=1:numel(stimID1);
        s = stimID1(i);
        idx = find(stm.outputsInfo.OutputStimID==s);
        o = stm.outputsInfo.OutputOrder(idx);
        if o==0
            powers(i) = 0;
        else
            powers(i) = powerList(o);
        end
    end
    
    %     window = [6:12];
    %
    %     trialsTouse = stm.lowMotionTrials & mean(stm.runVal,2)'<6;
    %
    %     datMean(:,1) = mean(mean(datAlign(:,window,trialsTouse & stimID==us(1)),3),2);
    %     for pidx = 1:numel(powerList);
    %         s = us(find(outputsOrder == pidx));
    %         datMean(:,pidx+1) =  mean(mean(datAlign(:,window,trialsTouse & stimID==s),3),2);
    %     end
    
    
    PSTHs = permute(datAlign,[3 1 2]);
    powers = powers;
    
    
    
    pbal = powerCalib(PSTHs,powers,0);
    pbal.calibCells = find(calibCells); 
    
    calibCellsList = find(calibCells);
    clim = [0 5];
    
    if plotExamples
        for i=1:numel(calibCellsList)
            if pbal.Rsquare(i) > 0.25 && pbal.cellsThatFit(i)
                    figure(1001)
                    clf;
                
                pl = unique(powers);
                nPwr = numel(pl);
                
                colorList = colorMapPicker(nPwr,'viridis');
                
                baseVal = datAlign(i,:,powers==0);
                baseVal = nanmean(baseVal(:)); 

                
                meanVal=[];seVal=[];
                for k=1:numel(pl)
                    ax(k) = subplot(2,nPwr+1,k);
                    cellDat = squeeze(datAlign(i,:,powers==pl(k))) - baseVal;
                    fillPlot(cellDat',[],'ci',colorList{k},'none',colorList{k},0.5);
                    box off
                    
                    subplot(2,nPwr+1,k+nPwr+1)
                    imagesc(cellDat')
                    caxis(clim)
                    
                    cellDat = cellDat(6:end,:);
                    meanVal(k) = nanmean(cellDat(:));
                    
                    trialDat = nanmean(cellDat,1);
                    seVal(k) = nanstd(trialDat)/sqrt(numel(trialDat));
                end
                subplot(2,nPwr+1,k+1)
                errorbar(pl,meanVal,seVal)
                hold on
                plot(pbal.allFits{i})
                xlim([0 .125])
                legend off
                title({['Ind: ' num2str(ind) ' cell: ' num2str(calibCellsList(i))],...
                    ['Power At Sat ' num2str(pbal.PowerAtSat(i))], ...
                    ['Fit R2: ' num2str(pbal.Rsquare(i))] } )
                linkaxes(ax)
                
                pause
            end
        end
    end
%     
    %     powerList
    %     satPowerForFitCells = cat(2,satPowerForFitCells,pbal.PowerAtSat(pbal.cellsThatFit));
    satPowerForFitCellsTemp{ind}= pbal.PowerAtSat(pbal.cellsThatFit);
    
    allPBal{ind}=pbal;
    %
    allPowerCurveData{ind}.PSTHs = PSTHs;
    allPowerCurveData{ind}.powers = powers;
end
satPowerForFitCells = cat(2,satPowerForFitCellsTemp{:});
toc
disp('done')
%% Some Plots for Power Calib
cellsThatFit=[];
for ind = 1:numExps;
    cellsThatFit = cat(2,cellsThatFit,allPBal{ind}.cellsThatFit);
end

R2=[];
for ind = 1:numExps;
    R2=cat(2,R2,allPBal{ind}.Rsquare(allPBal{ind}.cellsThatFit));
end

figure(93);clf
histogram(satPowerForFitCells,17)
ylabel('Count')
xlabel('Target Power (W)')

disp(['Mean : ' num2str(mean(satPowerForFitCells)*1000,3) ' ' char(177) ' ' num2str(ste(satPowerForFitCells)*1000,2) ' mW'])
disp(['N = ' num2str(numel(satPowerForFitCells))])


figure(94);clf;
scatter(satPowerForFitCells,R2)
ylabel('Fit R2')
xlabel('Target Power')

%% Plot Example Trace


% for ind =8
%        stm = All(ind).out.stm;
%     FR1 = All(ind).out.info.FR;
%     
%     try
%         stime = stm.holoRequest.bigListOfFirstStimTimes(:,1);
%     catch
%         stime = stm.holoRequest2.bigListOfFirstStimTimes(:,1);
%     end
%     stFr = round(stime*FR1);
%     
%     dataToUse1 = stm.zdfData;
%     sz = size(dataToUse1);
%     
%     trialsToUse = stm.lowMotionTrials & mean(stm.runVal,2)'<6;
%     dataToUse1(:,:,~trialsToUse)=nan;
%     
%     calibCells = allPBal{ind}.calibCells;
%     
% %     tc = stm.targetedCells;
% 
%     
% end

%% Spike Curve Plots

subtractNull=0;
baselineSubtract=1;

satPowerForFitCells =[];
allPBal =[];
counter2 = 0;

plotFast=1;
plotIt=0;0;

plotExample =0;

Rsquare2=[];
MDatMat=[];
coeff2 =[];

%Note only 3 to end were done the same way (correctly)
for ind = 1:numExps;
    spk = All(ind).out.spk;
    FR = All(ind).out.info.FR;
    
    stime = spk.holoRequest.bigListOfFirstStimTimes(:,1);
    stFr = round(stime*FR);
    
    dataToUse = spk.caiman_matched;    spk.zdfData;
    sz = size(dataToUse);
    
    
    trialsToUse = spk.lowMotionTrials & mean(spk.runVal,2)'<6;
    dataToUse(:,:,~trialsToUse)=nan;
    
    tc = spk.targetedCells;
    
    calibCells  = ~isnan(tc) & ~isnan(stime);
    afterFrames = 18;
    preFrames =5;
    
    datAlign = nan([sum(calibCells) preFrames+afterFrames+1 sz(3)]);
    counter =1;
    for i=1:numel(tc);
        if calibCells(i);
            c = tc(i);
            st = stFr(i);
            if st+afterFrames > sz(2)
                padSize = st+afterFrames-sz(2);
                
                datpart = dataToUse(c,st-preFrames:end,:);
                datAlign(counter,:,:) = padarray(datpart,[0 padSize 0],nan,'post');
                disp('Padding')
            else
                datAlign(counter,:,:) = dataToUse(c,st-preFrames:st+afterFrames,:);
            end
            counter=counter+1;
        end
    end
    
    powerList = spk.holoRequest.holoStimParams.powerList{:};
    spikeList = spk.holoRequest.holoStimParams.pulseList(:);
    
    try
        outputsOrder = spk.outputsInfo.OutputOrder;
    catch
        outputsOrder = spk.stimParams.Seq;
    end
    
    stimID = spk.stimID;
    %     stimID = stimID(1:77);
    us = unique(stimID);
    
    spikes=[];
    for i=1:numel(stimID);
        s = stimID(i);
        %         idx = find(spk.outputsInfo.OutputStimID==s);
        idx = find(us==s);
        
        o = spk.outputsInfo.OutputOrder(idx);
        if o==0
            spikes(i) = 0;
        else
            spikes(i) = spikeList(o);
        end
    end
    
    PSTHs = permute(datAlign,[3 1 2]);
    spikes = spikes;
    
    baselinePeriod =1:5;
    samplePeriod =9:18;
    
    sigThreshold = 2;
    
    spikeList = unique(spikes);
    
    c=0;
    for i=1:numel(spikeList)
        p=spikeList(i);
        numPassAvg = sum(spikes==p);
        x = PSTHs(spikes==p,:,:);
        m = squeeze(nanmean(x,1));
        
        x = permute(x,[2 1 3]);
        dataPeriods = nanmean(x(:,:,samplePeriod),3);
        dataPeriodsB = nanmean(x(:,:,baselinePeriod),3);
        
        if subtractNull
            NullData = PSTHs(spikes==0,:,:);
            NullData = permute(NullData,[2 1 3]);
            dataPeriodsNull = nanmean(NullData(:,:,samplePeriod),3);
            dataPeriods = bsxfun(@minus,dataPeriods,nanmean(dataPeriodsNull,2));
            dataPeriodsB = bsxfun(@minus,dataPeriodsB,nanmean(dataPeriodsNull,2));
        end
        if baselineSubtract
            dataPeriods = dataPeriods- dataPeriodsB;
        end
        mDat{i} = dataPeriods;
        [htst , pVals{i}] = ttest(dataPeriodsB', dataPeriods') ;% for stimmability
    end
    
    
    
    % coeff2 = [];
    % Rsquare2 = [];
    % fits=[];
    
    
    if plotIt
        figure(10);
    end
    
    % clear minStim maxStim minFluor maxFluor
    
    smallCounter=0;
    Rsquare2One=[];
    
    FitCells=1:size(mDat{1},1);
    for i=1:size(mDat{1},1);
        
        if mod(i,20)==0
            fprintf([num2str(i) ' ' ])
        end
        
        c=FitCells(i);
        counter2=counter2+1;
        smallCounter =smallCounter+1;
        if plotIt
            clf
        end
        tempDat = cellfun(@(x) x(c,:),mDat,'uniformoutput',0);
        mmDat = cellfun(@(x) nanmean(x(c,:)),mDat);
        sdDat = cellfun(@(x) nanstd(x(c,:)),mDat);
        nmDat = cellfun(@(x) numel(x(c,:)),mDat);
        seDat = sdDat./sqrt(nmDat);
        
        spksUsed = unique(spikes); %[ 0 1 3 10 30];
        tDat = [];
        gDat = [];
        for i=1:numel(tempDat);
            tDat = [tDat tempDat{i}];
            gDat = [gDat ones([1 numel(tempDat{i})])*spksUsed(i)];
        end
        
        
        if plotIt
            subplot(1,2,1)
            plotSpread(tempDat);
            e = errorbar(mmDat,seDat);
            e.LineWidth = 2;
            e.Color = rgb('grey');
            title(c);
            ylabel('Fluorescence')
            xlabel('Stim')
            
            
            subplot(1,2,2);
            scatter(gDat,tDat)
        end
        try
            excl = isnan(tDat) | isnan(gDat);
            gDat(excl)=[];
            tDat(excl)=[];
            ft =fittype({'x.^2','x','1'});
            
            
            [f2 gof2] = fit(gDat',tDat',ft);
            [f2inv gof2] = fit(tDat',gDat',ft);
            
            % y=tDat;
            % y = y-min(y);
            % x=gDat;
            % % x=x+1;
            %
            % subplot(1,2,1)
            % myfittype = fittype({'log(x+1)', '1'},'coefficients',{'a1','a2'})
            % [myfit gof] = fit(x',y',myfittype);
            % % myfit = fit(y',x',myfittype)
            %
            % plot(myfit,x,y,'o')
            
            % subplot(1,2,2)
            % myfittype = fittype({'exp(x+1)', '1'},'coefficients',{'a1','a2'})
            %  myfit = fit(y',x',myfittype)
            %
            % plot(myfit,y,x,'o')
            
            % plot(x,y,'o')
            
            if plotIt
                
                hold on
                e = errorbar(spksUsed,mmDat,seDat);
                e.LineWidth = 2;
                e.Color = rgb('grey');
                p2 = plot(f2);
                p2.Color ='k';
                p2.LineWidth=2;
            end
            
            coeff2(counter2) = f2inv.a; %f2.p1;
            Rsquare2(counter2) = gof2.rsquare;
            Rsquare2One(smallCounter) = gof2.rsquare;
            
            fits{counter2} = f2inv;
            fits2{counter2}= f2;
            
        catch
            coeff2(counter2) = nan;
            Rsquare2(counter2) =0;
            Rsquare2One(smallCounter) = 0;
            
            f = fittype('a*x^2 + b*x + c');
            fits{counter2} = cfit(f,0,0,0);
            %         disp('.')
        end
        if plotIt
            
            ylabel('Fluorescence')
            xlabel('Pulses added')
            
            zeroint=fits2{counter2}(0);
            title({['R2: ' num2str(Rsquare2(counter2))]; ['0 intercept: ' num2str(zeroint)]})
            
        end
        
        minStim(counter2) = mmDat(1);
        maxStim(counter2) = mmDat(end);
        
        if ~plotFast &  plotIt
            drawnow
            pause
        end
        
        if plotExample & Rsquare2(counter2)>0.75
        disp('plot?')
            figure(1);clf
            nspk = numel(spikeList);
            colorList = colorMapPicker(nspk+1,'magma');
            for k=1:numel(spikeList)
                p=spikeList(k);
                numPassAvg = sum(spikes==p);
                x = PSTHs(spikes==p,c,:); %order: trials cells frames
                x=squeeze(x); 
                
                ax(k) = subplot(1,nspk+1,k)
                fillPlot(x,[],'ci',colorList{k},'none',colorList{k},0.5);
                box off
                xticks(0:FR:24)
                xticklabels(0:4)
            end
            linkaxes(ax)
            xlim([0 24])
            ylim([-1 5])
            
            subplot(1,nspk+1,k+1);
            s = scatter(gDat,tDat);
            s.MarkerFaceColor = rgb('darkorange');
            s.MarkerEdgeColor = 'none';
            s.SizeData=20;
            
            
            ylabel('Fluorescence')
            xlabel('Pulses added')
            hold on
            e = errorbar(spksUsed,mmDat,seDat);
            e.LineWidth = 2;
            e.Color = rgb('grey');
            p2 = plot(f2);
            p2.Color ='k';
            p2.LineWidth=2;
            
            zeroint=fits2{counter2}(0);
            title({['ind: ' num2str(ind) ' cell ' num2str(c)], ['R2: ' num2str(Rsquare2(counter2))], ['0 intercept: ' num2str(zeroint)]})
            legend off
            
          pause
        end
    end
    
    mDatMat = cell2mat(cellfun(@(x) nanmean(x,2), mDat,'uniformOutput',0));
    MDatMat = cat(1,MDatMat,mDatMat); %grand all data set;
    
    includeThreshold = Rsquare2One>0.15;
    mDatMat(~includeThreshold,:)=[];
    sum(includeThreshold)
    
    figure(101);clf
    p = plot(spikeList,mDatMat,'color',rgb('grey'));
    hold on
    m = mean(mDatMat);
    sd = std(mDatMat);
    se = sd./sqrt(size(mDatMat,1));
    e = errorbar(spikeList,m,se);
    e.Color = rgb('black');
    e.LineWidth = 2;
    ylim([-1 5])
    
    if plotExample
        
        
    end
    
    
    
    pause
    disp('.')
    
    
    
end

% coef cellfun(@(x)x.a,fits2)
coeff2
sum(Rsquare2>0.15)

figure(103);clf
% scatter(Rsquare2,coeff2,'o')

includeThreshold = Rsquare2>0.15;
MDatMat(~includeThreshold,:)=[];
sum(includeThreshold)

p = plot(spikeList,MDatMat,'color',rgb('grey'));
hold on
m = mean(MDatMat);
sd = std(MDatMat);
se = sd./sqrt(size(MDatMat,1));
e = errorbar(spikeList,m,se);
e.Color = rgb('black');
e.LineWidth = 2;
xticks(0:5:20)


%% Establish PCA on Vis Epoch
numExps = numel(All);

stepByStep=0;
useBothVisAndHoloForPCA =0;
plotHoloMerge =1;
plotMerge = 1;
plotNoStim =0;

plotExamples=0;

pauseEachInd=1;

% ID which holo corresponds with which ID
hIDS{1} = [1 9 7 4 6];
hPlot{1} = ones([1 5]);

hIDS{2} = [1:9];
hPlot{2} =  [1 1 1 1 1 1 1 1 1 ];

hIDS{3} = [1 2 3 4 5];
hPlot{3} = ones([1 5]);

for i =4:numExps
    hIDS{i} = [1 1:9];
    hPlot{i} = [1 0 1 1 1 1 1 1 1 1 ];
end


cosSimOut=[];
cosSimOutShuffle=[];
cosSimOnDiagOut=[];
cosSimOnDiagShuffleOut=[];
cosSimOffDiagOut=[];
CosSimToSpont=[];

rangeToAnalyze =1:24; %for PCA
plotRange = 1:15;
smoothFactor =3;

   numCellsPerManifold = [];
    numCellsPerPattern = [];
    numSpikesPerPattern = [];
    
allVisTC=[];
allHoloTC=[];

for ind = 1:numExps;
    visID = All(ind).out.vis.visID;
    uv = unique(All(ind).out.vis.visID);
    
    try
        rtgs = unique(cat(2,All(ind).out.mani.rois{:}));
    catch
        rtgs = unique(cat(1,All(ind).out.mani.rois{:}));
    end
    htgs = All(ind).out.exp.targetedCells(rtgs);
    htgs = rtgs;
    htgs(isnan(htgs))=[];
    cellsToUse = htgs;
    
    %data to build PCA
    %build new composite zdfData;
    dat1 =All(ind).out.vis.caiman; %All(ind).out.vis.caiman_matched; All(ind).out.vis.allData;
    dat2 = All(ind).out.mani.caiman;%All(ind).out.mani.caiman_matched; All(ind).out.mani.allData;
    
    sz1 = size(dat1);
    sz2 = size(dat2);
    mnFrames = min(sz1(2),sz2(2));
    
    
    fullDat = cat(3,dat1(:,1:mnFrames,:),dat2(:,1:mnFrames,:));
    [fulldfDat, fullzdfDat] = computeDFFwithMovingBaseline(fullDat);
    
    dat1 = fullzdfDat(:,:,1:sz1(3));
    dat2 = fullzdfDat(:,:,sz1(3)+1:end);
    
    
    
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
            sz = size(datToPlot);
            %    datToPlot = smooth(datToPlot,3);
            datToPlot = reshape(datToPlot,sz);
            %         p = plot3(datToPlot(:,1),datToPlot(:,2),datToPlot(:,3));
            p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
            p.Color = colors{i};
            p.LineWidth = 2;
            
        else
            if i>1 && i<6
                datToPlot = squeeze(mean(visScore(:,plotRange,(visID==uv(i) | visID==uv(i+4)) & trialsToUseVis),3));
                p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                p.Color = colors{i};
                p.LineWidth = 2;
            end
        end
        
        if stepByStep
            pause
        end
        
    end
    %     clf;hold on
    %     us = All(ind).out.mani.uniqueStims;
    %     stimID = All(ind).out.mani.stimID;
    us = All(ind).out.mani.uniqueStims;
    
    stimID = All(ind).out.mani.stimID;
    [s sidx] = sort(All(ind).out.mani.stimParams.Seq);
    us = us(sidx);
    
    if numel(us)<9
        us = [us 0 0 0 0 0 0 0 0 0]; %trailing zeros shouldn't matter
    end
    
    
    holoIDS = [1 1:9];
    holoPlot = [1 0 1 1 1 1 1 1 1 1 ];
    
    %     holoIDS = hIDS{ind};
    %     holoPlot = hPlot{ind};
    
    
    %     datToPlot = squeeze(mean(newHScore(:,:,stimID==us(2) & trialsToUseExp),3));
    %     p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
    %     p.Color = colors{1}; %rgb('red'); %colors{i};
    %     p.LineWidth =2;
    %     p.LineStyle = ':';
    
    
    for i=1:numel(us)
        if us(i)~=0
            if ~plotHoloMerge
                if holoIDS(i)==1 && ~plotNoStim
                    %do nothing
                else
                    datToPlot = squeeze(mean(newHScore(:,plotRange,stimID==us(i) & trialsToUseExp),3));
                    p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                    p.Color = colors{holoIDS(i)}; %rgb('red'); %colors{i};
                    p.LineStyle = ':';
                    if holoIDS(i)>=6
                        p.LineWidth =3;
                    else
                        p.LineWidth =2;
                    end
                end
            else
                if ~holoPlot(i)
                    %do nothing
                elseif holoIDS(i)>1 & holoIDS(i)<6; %i>2 && i<7
                    newIDS = holoIDS(i)+4;
                    tempID = us(find(holoIDS==newIDS));
                    if isempty(tempID)
                        tempID=NaN;
                    end
                    datToPlot = squeeze(mean(newHScore(:,plotRange,(stimID==us(i) | stimID==tempID)  & trialsToUseExp),3));
                    p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                    p.Color = colors{holoIDS(i)}; %rgb('red'); %colors{i};
                    p.LineWidth =2;
                    p.LineStyle = ':';
                elseif holoIDS(i)>=6; %i>2 && i<7
                    newIDS = holoIDS(i)-4;
                    tempID = us(find(holoIDS==newIDS));
                    if isempty(tempID)
                        tempID=NaN;
                    end
                    datToPlot = squeeze(mean(newHScore(:,plotRange,(stimID==us(i) | stimID==tempID)  & trialsToUseExp),3));
                    p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                    p.Color = colors{holoIDS(i)}; %rgb('red'); %colors{i};
                    p.LineWidth =2;
                    p.LineStyle = ':';
                elseif holoIDS(i)==1 && plotNoStim
                    datToPlot = squeeze(mean(newHScore(:,plotRange,stimID==us(i) & trialsToUseExp),3));
                    p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                    p.Color = colors{holoIDS(i)}; %rgb('red'); %colors{i};
                    p.LineWidth =2;
                    p.LineStyle = ':';
                end
            end
            
            if stepByStep
                pause
            end
        end
    end
    
    title(['Ind: ' num2str(ind)]);
    xlabel(['PC' num2str(ax1)])
    ylabel(['PC' num2str(ax2)])
    zlabel(['PC' num2str(ax3)])
    %     pause
    
    %% Cosine Similarity by Orientation
    figure(4);clf
    subplot(1,3,1)
    colormap parula
    
    crange = [0 1];
    respRange =6:18;
    
    %
    %     dat1 = All(ind).out.vis.zdfData;
    %     dat2 = All(ind).out.exp.zdfData;
    %
    %     sz1 = size(dat1);
    %     sz2 = size(dat2);
    %     mnFrames = min(sz1(2),sz2(2));
    %
    %
    %     fullDat = cat(3,All(ind).out.vis.allData(:,1:mnFrames,:),All(ind).out.exp.allData(:,1:mnFrames,:));
    %     [fulldfDat fullzdfDat] = computeDFFwithMovingBaseline(fullDat);
    % %
    dat1 = fullzdfDat(:,:,1:sz1(3));
    dat2 = fullzdfDat(:,:,sz1(3)+1:end);
    %     dataToUse = All(ind).out.vis.zdfData(cellsToUse,:,:);
    %     hDataToUse = All(ind).out.exp.zdfData(cellsToUse,:,:);
    
    
    dataToUse = dat1(cellsToUse,:,:);
    hDataToUse = dat2(cellsToUse,:,:);
    
    bhData = mean(hDataToUse(:,1:5,:),2);
    hDataToUse = hDataToUse-bhData; %baseline
    %     hDataToUse = smoothdata(hDataToUse,2); %smooth
    
    bData = nanmean(dataToUse(:,1:5,:),2);
    %     dataToUse = dataToUse-bData; %baseline
    
    temp = permute(dataToUse,[1 3 2]);
    temp2 = squeeze(bData);
    temp3 = temp-temp2;
    dataToUse = permute(temp3,[1 3 2]);
    
    
    %     dataToUse = smoothdata(dataToUse,2); %smooth
    
    %     visMean-0.2zeros([numel(cellsToUse) 9]);
    %     for i=1:9
    %         val = mean(mean(datToUse(:,6:15,visID==uv(i)& trialsToUseVis),2),3);
    %         visMean(:,i)= val;
    %     end
    
    visMean=zeros([numel(cellsToUse) 5]);
    val = mean(mean(dataToUse(:,respRange,visID==uv(1)& trialsToUseVis),2),3);
    visMean(:,1)= val;
    for i=2:5
        val = mean(mean(dataToUse(:,respRange,(visID==uv(i)| visID==uv(i+4)) & trialsToUseVis),2),3);
        visMean(:,i)= val;
    end
    
    
    imagesc(visMean)
    title('Vis Mean Response')
    colorbar
    caxis(crange)
    
    subplot(1,3,2)
    
    %     holoMean = zeros([numel(cellsToUse) numel(us)]);
    %     for i=1:numel(us)
    %         val = mean(mean(hdatToUse(:,6:15,stimID==us(i)& trialsToUseExp),2),3);
    %         holoMean(:,i)= val;
    %     end
    holoMean = zeros([numel(cellsToUse) 5]);
    val = mean(mean(hDataToUse(:,respRange,stimID==us(1)& trialsToUseExp),2),3);
    holoMean(:,1)= val;
    %     val = mean(mean(hDataToUse(:,respRange,stimID==us(2)& trialsToUseExp),2),3);
    %     holoMean(:,2)= val;
    for i= find(holoIDS>1 & holoIDS<=5); %3:6
        hi = holoIDS(i);
        newIDS = holoIDS(i)+4;
        tempID = us(find(holoIDS==newIDS));
        if isempty(tempID)
            tempID=NaN;
        end
        val = mean(mean(hDataToUse(:,respRange,(stimID==us(i) | stimID==tempID) & trialsToUseExp),2),3);
        holoMean(:,hi)= val;
    end
    imagesc(holoMean)
    title('Holo Mean Response')
    colorbar
    caxis(crange)
    
    cosSim=zeros([size(visMean,2) size(holoMean,2)]);
    for i=1:size(visMean,2)
        for k=1:size(holoMean,2)
            cosSim(i,k)=cosine_similarity(visMean(:,i),holoMean(:,k));
        end
    end
    subplot(1,3,3)
    imagesc(cosSim)
    colorbar
    
    %
    %     onDiag =[12 18 24 30];% [6 12 18 24 30];
    %     meanSimOnDiag = mean(cosSim(onDiag));
    %     onDiag = [onDiag 1:11 16 21 26];
    %     meanSimOffDiag = mean(cosSim(~ismember(1:numel(cosSim),onDiag)));
    %     meanSimOnDiag/meanSimOffDiag
    
    
    testSim = cosSim(2:5,2:5);
    
    meanSimOnDiag = mean(testSim(logical(eye(4))));
    meanSimOffDiag = mean(testSim(~logical(eye(4))));
    disp(['Holography (ori) mean on diagonal : ' num2str(meanSimOnDiag) ', off diag: ' num2str(meanSimOffDiag)])
    
    title('Cosine Similarity, Orientation')
    
    %% cosine similarity by direction not ori
    
    figure(5);clf
    colormap parula
    
    crange = [0 2];
    respRange =6:18;%9:18;
    
    subplot(1,3,1)
    visMean=zeros([numel(cellsToUse) 9]);
    val = mean(mean(dataToUse(:,respRange,visID==uv(1)& trialsToUseVis),2),3);
    visMean(:,1)= val;
    for i=2:9
        val = mean(mean(dataToUse(:,respRange,(visID==uv(i)) & trialsToUseVis),2),3);
        visMean(:,i)= val;
    end
    
    
    imagesc(visMean)
    title('Vis Mean Response')
    colorbar
    caxis(crange)
    
    subplot(1,3,2)
    
    %     holoMean = zeros([numel(cellsToUse) numel(us)]);
    %     for i=1:numel(us)
    %         val = mean(mean(hdatToUse(:,6:15,stimID==us(i)& trialsToUseExp),2),3);
    %         holoMean(:,i)= val;
    %     end
    holoMean = zeros([numel(cellsToUse) 9]);
    val = mean(mean(hDataToUse(:,respRange,stimID==us(1)& trialsToUseExp),2),3);
    holoMean(:,1)= val;
    %     val = mean(mean(hDataToUse(:,respRange,stimID==us(2)& trialsToUseExp),2),3);
    %     holoMean(:,2)= val;
    for i= find(holoIDS>1); %3:6
        hi = holoIDS(i);
        val = mean(mean(hDataToUse(:,respRange,(stimID==us(i) ) & trialsToUseExp),2),3);
        holoMean(:,hi)= val;
    end
    imagesc(holoMean)
    title('Holo Mean Response')
    colorbar
    caxis(crange)
    
    cosSim=zeros([size(visMean,2) size(holoMean,2)]);
    for i=1:size(visMean,2)
        for k=1:size(holoMean,2)
            cosSim(i,k)=cosine_similarity(visMean(:,i),holoMean(:,k));
        end
    end
    subplot(1,3,3)
    imagesc(cosSim)
    colorbar
    
    
    testSim = cosSim(1:9,2:9);
    cosSimOut = cat(2,cosSimOut,testSim);
    size(cosSimOut);
    
    testSim = cosSim(2:9,2:9);
    meanSimOnDiag = nanmean(testSim(logical(eye(8))));
    meanSimOffDiag = nanmean(testSim(~logical(eye(8))));
    meanSimOnDiag/meanSimOffDiag;
    cosSimOnDiagOut = cat(1,cosSimOnDiagOut,testSim(logical(eye(8))));
    cosSimOffDiagOut = cat(1,cosSimOffDiagOut,testSim(~logical(eye(8))));
    CosSimToSpont  = cat(2,CosSimToSpont, cosSim(1,2:9));
    
    
    disp(['Holography mean on diagonal : ' num2str(meanSimOnDiag) ', off diag: ' num2str(meanSimOffDiag)])
    
    yticks(1:9)
    yticklabels({'gray' 0:45:315})
    ylabel('Visual Input')
    xticks(1:10)
    xticklabels({'no stim' 0:45:315})
    xtickangle(45)
    xlabel('Holographic Input')
    %     title('Cosine Similarity')
    title('Cosine Similarity, Direction')
    
    caxis([-0.2 0.75])
    
    %% shuffle Cosine Similarity
    
    figure(15);clf
    colormap parula
    
    crange = [0 2];
    respRange =6:18;%9:18;
    
    shuffleOrder = randperm(size(dataToUse,3));
    dataToUseShuffle = dataToUse(:,:,shuffleOrder);
    
    subplot(1,3,1)
    visMean=zeros([numel(cellsToUse) 9]);
    val = mean(mean(dataToUseShuffle(:,respRange,visID==uv(1)& trialsToUseVis),2),3);
    visMean(:,1)= val;
    for i=2:9
        val = mean(mean(dataToUseShuffle(:,respRange,(visID==uv(i)) & trialsToUseVis),2),3);
        visMean(:,i)= val;
    end
    
    
    imagesc(visMean)
    title('Vis Mean Response')
    colorbar
    caxis(crange)
    
    subplot(1,3,2)
    
    %     holoMean = zeros([numel(cellsToUse) numel(us)]);
    %     for i=1:numel(us)
    %         val = mean(mean(hdatToUse(:,6:15,stimID==us(i)& trialsToUseExp),2),3);
    %         holoMean(:,i)= val;
    %     end
    shuffleOrder = randperm(size(hDataToUse,3));
    hDataToUseShuffle =  hDataToUse(:,:,shuffleOrder);
    
    holoMean = zeros([numel(cellsToUse) 9]);
    val = mean(mean(hDataToUseShuffle(:,respRange,stimID==us(1)& trialsToUseExp),2),3);
    holoMean(:,1)= val;
    %     val = mean(mean(hDataToUse(:,respRange,stimID==us(2)& trialsToUseExp),2),3);
    %     holoMean(:,2)= val;
    for i= find(holoIDS>1); %3:6
        hi = holoIDS(i);
        val = mean(mean(hDataToUseShuffle(:,respRange,(stimID==us(i) ) & trialsToUseExp),2),3);
        holoMean(:,hi)= val;
    end
    imagesc(holoMean)
    title('Holo Mean Response')
    colorbar
    caxis(crange)
    
    cosSim=zeros([size(visMean,2) size(holoMean,2)]);
    for i=1:size(visMean,2)
        for k=1:size(holoMean,2)
            cosSim(i,k)=cosine_similarity(visMean(:,i),holoMean(:,k));
        end
    end
    subplot(1,3,3)
    imagesc(cosSim)
    colorbar
    
    
    testSim = cosSim(1:9,2:9);
    cosSimOutShuffle = cat(2,cosSimOutShuffle,testSim);
    
    testSim = cosSim(2:9,2:9);
    meanSimOnDiag = nanmean(testSim(logical(eye(8))));
    meanSimOffDiag = nanmean(testSim(~logical(eye(8))));
    meanSimOnDiag/meanSimOffDiag;
    cosSimOnDiagShuffleOut = cat(1,cosSimOnDiagShuffleOut,testSim(logical(eye(8))));
    
    disp(['Holography Shuffled mean on diagonal : ' num2str(meanSimOnDiag) ', off diag: ' num2str(meanSimOffDiag)])
    
    yticks(1:9)
    yticklabels({'gray' 0:45:315})
    ylabel('Visual Input')
    xticks(1:10)
    xticklabels({'no stim' 0:45:315})
    xtickangle(45)
    xlabel('Holographic Input')
    %     title('Cosine Similarity')
    title('Cosine Similarity, Direction')
    
    caxis([-0.2 0.75])
    
    %% what does cosine similarity look like within vis trials
    crange = [0 2];
    %     respRange =6:18;
    
    halfTrials = zeros(size(trialsToUseVis));
    %     oddTrial(1:2:end) = 1;
    halfTrials(randperm(numel(halfTrials),floor(numel(halfTrials)/2))) =1;
    
    
    figure(6);clf
    subplot(1,3,1)
    visMean1=zeros([numel(cellsToUse) 9]);
    val = mean(mean(dataToUse(:,respRange,visID==uv(1)& trialsToUseVis & ~halfTrials ),2),3);
    visMean1(:,1)= val;
    for i=2:9
        val = mean(mean(dataToUse(:,respRange,(visID==uv(i)) & trialsToUseVis & ~halfTrials ),2),3);
        visMean1(:,i)= val;
    end
    imagesc(visMean1)
    title('Vis Mean Response')
    colorbar
    caxis(crange)
    
    subplot(1,3,2)
    visMean2=zeros([numel(cellsToUse) 9]);
    val = mean(mean(dataToUse(:,respRange,visID==uv(1)& trialsToUseVis & halfTrials),2),3);
    visMean1(:,1)= val;
    for i=2:9
        val = mean(mean(dataToUse(:,respRange,(visID==uv(i)) & trialsToUseVis & halfTrials),2),3);
        visMean2(:,i)= val;
    end
    imagesc(visMean2)
    title('Vis Mean Response')
    colorbar
    caxis(crange)
    
    
    cosSim=zeros([size(visMean1,2) size(visMean2,2)]);
    for i=1:size(visMean1,2)
        for k=1:size(visMean2,2)
            cosSim(i,k)=cosine_similarity(visMean1(:,i),visMean2(:,k));
        end
    end
    subplot(1,3,3)
    imagesc(cosSim)
    colorbar
    
    
    testSim = cosSim(2:9,2:9);
    
    meanSimOnDiag = nanmean(testSim(logical(eye(8))));
    meanSimOffDiag = nanmean(testSim(~logical(eye(8))));
    meanSimOnDiag/meanSimOffDiag;
    disp(['Vis mean on diagonal : ' num2str(meanSimOnDiag) ', off diag: ' num2str(meanSimOffDiag)])
    
    
    yticks(1:9)
    yticklabels({'gray' 0:45:315})
    ylabel('Visual Input')
    xticks(1:9)
    xticklabels({'gray' 0:45:315})
    xtickangle(45)
    xlabel('Visual Input')
    title('Cosine Similarity')
    ylabel('Visual Input')
    title('Cosine Similarity, Visual')
    
    caxis([-0.2 0.75])
    
    %% figure out Cells and spikes per stim
    
    cellsPerROI = cellfun(@numel,All(ind).out.mani.rois);
    
    allcellTargets = [];
    
    for i=2:numel(us);
        r = find(All(ind).out.mani.stimParams.Seq == i-1) ;
        
        roiPointer = All(ind).out.mani.stimParams.roi{r};
        roiPointer = unique([All(ind).out.mani.holoRequest.bigListofSequences{i-1}]);
        if roiPointer==0
            cellsHitPerPattern(i)=0;
        else
            try
                cellsTargeted = unique(cat(1,All(ind).out.mani.rois{roiPointer}));
            catch
                cellsTargeted = unique(cat(2,All(ind).out.mani.rois{roiPointer}))';
            end
            cellsHitPerPattern(i) = numel(cellsTargeted);
            allcellTargets = [allcellTargets cellsTargeted'];
            
            stimPattern = ([All(ind).out.mani.holoRequest.bigListofSequences{i-1}]);
            totalAPsPerPattern(i) = sum(cellsPerROI(stimPattern));
        end
    end
    
    numCellsPerManifold(ind) = numel(unique(allcellTargets));
    numCellsPerPattern(ind,:) = cellsHitPerPattern;
    numSpikesPerPattern(ind,:) = totalAPsPerPattern;
    
    
    %% calc individual cell similarity
    visTC=[];holoTC=[];
    for i =1:numel(cellsToUse)
          evalRange = 9:18;
            
            vMean=[];vSTD=[];vSEM=[];vN=[];
            for k=1:9
                tempDat = squeeze(mean(dataToUse((i),evalRange,visID==uv(k) & trialsToUseVis),2));
                sz = size(tempDat);
                vMean(k)=nanmean(tempDat);
                vSTD(k) = nanstd(tempDat);
                vN(k) = numel(tempDat);
                
            end
            vSEM = vSTD./(sqrt(vN));
            
            hMean=[];hSTD=[];hSEM=[];hN=[];
            for k=1:numel(us)
                tempDat = squeeze(mean(hDataToUse((i),evalRange,stimID==us(k) & trialsToUseExp),2));
                sz = size(tempDat);
                hMean(k)=nanmean(tempDat);
                hSTD(k) = nanstd(tempDat);
                hN(k) = numel(tempDat);
                
            end
            hSEM = hSTD./(sqrt(hN));
        
            
            visTC(i,:) = vMean;
            holoTC(i,:) = hMean([1 3:10]);
            
    end
    
    allVisTC{ind}=visTC;
    allHoloTC{ind} = holoTC; 
    
    %% Plot Example Cell
    if plotExamples
        for i =1:numel(cellsToUse)
            evalRange = 9:18;
            
            vMean=[];vSTD=[];vSEM=[];vN=[];
            for k=1:9
                tempDat = squeeze(mean(dataToUse((i),evalRange,visID==uv(k) & trialsToUseVis),2));
                sz = size(tempDat);
                vMean(k)=nanmean(tempDat);
                vSTD(k) = nanstd(tempDat);
                vN(k) = numel(tempDat);
                
            end
            vSEM = vSTD./(sqrt(vN));
            
            hMean=[];hSTD=[];hSEM=[];hN=[];
            for k=1:numel(us)
                tempDat = squeeze(mean(hDataToUse((i),evalRange,stimID==us(k) & trialsToUseExp),2));
                sz = size(tempDat);
                hMean(k)=nanmean(tempDat);
                hSTD(k) = nanstd(tempDat);
                hN(k) = numel(tempDat);
                
            end
            hSEM = hSTD./(sqrt(hN));
            
            
            figure(7);clf
            subplot(1,2,1)
            e2 = errorbar(1:9,vMean,vSEM);
            xticks(1:9)
            xticklabels({'gray' 0:45:315})
            
            hold on;
            e2 = errorbar(1:9,hMean([1 3:10]),hSEM([1 3:10]));
            %          xticks(1:9)
            %          xticklabels({'gray' 0:45:315})
            ylabel('zdF/F')
            xlabel('Stimulus')
            
            subplot(1,2,2)
            
            minVal = min([hMean vMean]);
            if minVal>0
                minVal=0;
            end
            
            rads = 0:pi/4:2*pi;
            polarplot(rads,[vMean(2:9) vMean(2)]-minVal);
            hold on
            polarplot(rads,[hMean(3:10) hMean(3)]-minVal);
            
            
            colors = colorMapPicker(5,'plasma');
            colors = colors(1:4);
            colors = cat(2,{rgb('grey')},colors,colors);
            
            try
                estSpikes = All(ind).out.mani.mani.estSpikes;
            catch
                estSpikes = All(ind).out.mani.estSpikes;
            end
            
            cellsToWrite = All(ind).out.mani.mani.CellsToWrite; %defines estSpikes in units of fitCells
            cellsToFit = All(ind).out.mani.mani.CellsToFit;
            
            cellIDs = cellsToFit(cellsToWrite);
            
            tc = All(ind).out.mani.targetedCells;
            cellIDsS2P = tc(cellIDs);
            
            idxInEstSpikes = find(cellIDsS2P==cellsToUse(i));
            
            try
                spikesThisCell = estSpikes(idxInEstSpikes,:);
            catch
                spikesThisCell=NaN([1 9]);
            end
            %         tcsofCellIds = All(ind).out.mani.targetedCells(All(ind).out.mani.CellIDs);
            %         tempi = find(tcsofCellIds==cellsToUse(i));
            %         writeCell = All(ind).out.mani.CellIDs(tempi);
            %     fitID = find(All(ind).out.mani.mani.CellsToFit==writeCell);
            %
            %
            %     idxInEstSpikes = find(estSpikeCells==cellsToUse(i));
            %     spikesThisCell = estSpikes(idxInEstSpikes,:);
            
            figure(8);clf
            
            ax=[];
            ax2=[];
            for k=1:9
                ax(end+1) = subplot(4,9,k+9);
                tempDat = squeeze(dataToUse(i,:,visID==uv(k)));
                imagesc(tempDat')
                caxis([-1 4])
                xticks(0:FR:24)
                xticklabels(-1:1:3)
                
                ax2(end+1) = subplot(4,9,k);
                fillPlot(tempDat',[],'ci',colors{k},[],colors{k},0.5);
                xticks(0:FR:24)
                xticklabels(-1:1:3)
                box off
            end
            
            holoToPlotIDX = [1 3:10];
            for k=1:numel(holoToPlotIDX)
                ax(end+1) = subplot(4,9,k+27);
                tempDat = squeeze(hDataToUse(i,:, stimID==us(holoToPlotIDX(k)) & trialsToUseExp));
                imagesc(tempDat')
                caxis([-1 4])
                xticks(0:FR:24)
                xticklabels(-1:1:3)
                
                ax2(end+1) = subplot(4,9,k+18);
                fillPlot(tempDat',[],'ci',colors{k},[],colors{k},0.5);
                xticks(0:FR:24)
                xticklabels(-1:1:3)
                title(num2str(spikesThisCell(k)))
                box off
            end
            linkaxes(ax2)
            ylim([-1 3])
            xlim([1 24])
            same_color_scale(ax)
            colormap viridis
            
            pause
            
        end
    end
    
    if pauseEachInd
        pause
    end
    
end

pOnvShuffle = ranksum(cosSimOnDiagOut,cosSimOnDiagShuffleOut);
pOnvOff = ranksum(cosSimOnDiagOut,cosSimOffDiagOut);
pOnvSpont = ranksum(cosSimOnDiagOut,CosSimToSpont);

disp(['P on diag vs off: ' num2str(pOnvOff)])
disp(['P on diag vs Spont: ' num2str(pOnvSpont)])
disp(['P on diag vs Shuffle: ' num2str(pOnvShuffle)])

figure(17);clf
plotSpread({cosSimOnDiagOut cosSimOffDiagOut CosSimToSpont cosSimOnDiagShuffleOut })


%% Plot indiv cell Tuning Curve Similarities

visTCs = cat(1,allVisTC{:});
holoTCs = cat(1,allHoloTC{:});

TCsim=[];
for i =1:size(visTCs,1)
    TCsim(i) = cosine_similarity(visTCs(i,:)',holoTCs(i,:)');
    
    TCsim(i) = circ_var(abs(visTCs(i,2:end)-holoTCs(i,2:end))');
end

TCsimShuffle=[];
for i =1:size(visTCs,1)
    tempVis = visTCs(i,:);
    tempVis = tempVis(randperm(9));
    TCsimShuffle(i) = cosine_similarity(tempVis',holoTCs(i,:)');
    
TCsimShuffle(i) = cosine_similarity(visTCs(randperm(size(visTCs,1),1),:)',holoTCs(i,:)');

TCsimShuffle(i) = circ_var(abs(visTCs(randperm(size(visTCs,1),1),2:end)- holoTCs(i,2:end))');

end

figure(18);clf
plotSpread({TCsim TCsimShuffle},'showMM',4)
ranksum(TCsimShuffle,TCsim)
%% Plot Data

cosSimToPlot=zeros(size(cosSimOut'));
for i =1:size(cosSimOut,2)
    thisCosSim = cosSimOut(:,i);
    [mx p] = max(thisCosSim(2:end));
    cosSimToPlot(i,:) = [thisCosSim(1) circshift(thisCosSim(2:end),-1*(p-1))'];
end

figure(8);clf
subplot(2,1,1)
imagesc(cosSimToPlot)
caxis([-0.2 0.75])
colorbar
xticks(1:9)
xticklabels({'Grey' 0:45:315})

subplot(2,1,2)

plot(mean(cosSimToPlot))


%% specficically Ind 5 code

ind = 5;
figure(5);clf
colormap parula

crange = [0 2];
respRange =6:18;%9:18;

subplot(1,3,1)
visMean=zeros([numel(cellsToUse) 9]);
val = mean(mean(dataToUse(:,respRange,visID==uv(1)& trialsToUseVis),2),3);
visMean(:,1)= val;
for i=2:9
    val = mean(mean(dataToUse(:,respRange,(visID==uv(i)) & trialsToUseVis),2),3);
    visMean(:,i)= val;
end


imagesc(visMean)
title('Vis Mean Response')
colorbar
caxis(crange)

subplot(1,3,2)

%     holoMean = zeros([numel(cellsToUse) numel(us)]);
%     for i=1:numel(us)
%         val = mean(mean(hdatToUse(:,6:15,stimID==us(i)& trialsToUseExp),2),3);
%         holoMean(:,i)= val;
%     end
holoMean = zeros([numel(cellsToUse) 10]);
val = mean(mean(hDataToUse(:,respRange,stimID==us(1)& trialsToUseExp),2),3);
holoMean(:,1)= val;
%     val = mean(mean(hDataToUse(:,respRange,stimID==us(2)& trialsToUseExp),2),3);
%     holoMean(:,2)= val;
for i=3:10
    val = mean(mean(hDataToUse(:,respRange,(stimID==us(i) ) & trialsToUseExp),2),3);
    holoMean(:,i)= val;
end
imagesc(holoMean)
title('Holo Mean Response')
colorbar
caxis(crange)

cosSim=zeros([size(visMean,2) size(holoMean,2)]);
for i=1:size(visMean,2)
    for k=1:size(holoMean,2)
        cosSim(i,k)=cosine_similarity(visMean(:,i),holoMean(:,k));
    end
end
subplot(1,3,3)
imagesc(cosSim)
colorbar


testSim = cosSim(2:9,3:10);

meanSimOnDiag = mean(testSim(logical(eye(8))));
meanSimOffDiag = mean(testSim(~logical(eye(8))));
meanSimOnDiag/meanSimOffDiag;

disp(['Holography mean on diagonal : ' num2str(meanSimOnDiag) ', off diag: ' num2str(meanSimOffDiag)])

yticks(1:9)
yticklabels({'gray' 0:45:315})
ylabel('Visual Input')
xticks(1:10)
xticklabels({'no stim' 'spont' 0:45:315})
xtickangle(45)
xlabel('Holographic Input')
%     title('Cosine Similarity')
title('Cosine Similarity, Direction')

caxis([-0.2 0.75])

%%what does cosine similarity look like within vis trials
crange = [0 2];
%     respRange =6:18;

halfTrials = zeros(size(trialsToUseVis));
%     oddTrial(1:2:end) = 1;
halfTrials(randperm(numel(halfTrials),floor(numel(halfTrials)/2))) =1;


figure(6);clf
subplot(1,3,1)
visMean1=zeros([numel(cellsToUse) 9]);
val = mean(mean(dataToUse(:,respRange,visID==uv(1)& trialsToUseVis & ~halfTrials ),2),3);
visMean1(:,1)= val;
for i=2:9
    val = mean(mean(dataToUse(:,respRange,(visID==uv(i)) & trialsToUseVis & ~halfTrials ),2),3);
    visMean1(:,i)= val;
end
imagesc(visMean1)
title('Vis Mean Response')
colorbar
caxis(crange)

subplot(1,3,2)
visMean2=zeros([numel(cellsToUse) 9]);
val = mean(mean(dataToUse(:,respRange,visID==uv(1)& trialsToUseVis & halfTrials),2),3);
visMean1(:,1)= val;
for i=2:9
    val = mean(mean(dataToUse(:,respRange,(visID==uv(i)) & trialsToUseVis & halfTrials),2),3);
    visMean2(:,i)= val;
end
imagesc(visMean2)
title('Vis Mean Response')
colorbar
caxis(crange)


cosSim=zeros([size(visMean1,2) size(visMean2,2)]);
for i=1:size(visMean1,2)
    for k=1:size(visMean2,2)
        cosSim(i,k)=cosine_similarity(visMean1(:,i),visMean2(:,k));
    end
end
subplot(1,3,3)
imagesc(cosSim)
colorbar


testSim = cosSim(2:9,2:9);

meanSimOnDiag = nanmean(testSim(logical(eye(8))));
meanSimOffDiag = nanmean(testSim(~logical(eye(8))));
meanSimOnDiag/meanSimOffDiag;
disp(['Vis mean on diagonal : ' num2str(meanSimOnDiag) ', off diag: ' num2str(meanSimOffDiag)])


yticks(1:9)
yticklabels({'gray' 0:45:315})
ylabel('Visual Input')
xticks(1:9)
xticklabels({'gray' 0:45:315})
xtickangle(45)
xlabel('Visual Input')
title('Cosine Similarity')
ylabel('Visual Input')
title('Cosine Similarity, Visual')

caxis([-0.2 0.75])