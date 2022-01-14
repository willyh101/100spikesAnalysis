clear
close all

masterTic = tic;

addpath(genpath('100spikesAnalysis'))
%% loadLists

ai203ManifoldLoadList;
% % allLoadList;

% loadPath = 'path/to/outfiles/directory';
loadPath = 'T:\Outfiles';

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

   All(ind).out.exp.dataToUse = All(ind).out.exp.zdfData;
   All(ind).out.vis.dataToUse = All(ind).out.vis.zdfData;
   
end
%%
opts.visRecWinRange = [0.5 1.5]; [0.5 1.5];
opts.runThreshold = 6;
opts.runValPercent = 0.75; %percent of frames that need to be below run threshold

All = cleanDataVisOnly(All,opts);

%%
disp('Calculating Vis stuff...')
opts.visAlpha = 0.05;
[All, outVars] = CalcPVisRFromVis(All,opts);
visPercent = outVars.visPercent;
outVars.visPercentFromVis = visPercent;

disp('ready')


%% Establish PCA on Vis Epoch
numExps = numel(All);

stepByStep=0;
useBothVisAndHoloForPCA =1;

plotExamples=0;


% ID which holo corresponds with which ID
 hIDS{1} = [1 9 7 4 6];
 hPlot{1} = ones([1 5]); 

 hIDS{2} = [1:9];
 hPlot{2} =  [1 1 1 1 1 1 1 1 1 ];

  hIDS{3} = [1 2 3 4 5];
 hPlot{3} = ones([1 5]); 
 
    hIDS{4} = [1 1:9];
    hPlot{4} = [1 0 1 1 1 1 1 1 1 1 ];
    
    hIDS{5} = [1 1:9];
    hPlot{5} = [1 0 1 1 1 1 1 1 1 1 1 ];

    cosSimOut=[];

for ind = 4:5; 1:numExps;
    visID = All(ind).out.vis.visID;
    uv = unique(All(ind).out.vis.visID);
    
    try
        rtgs = unique(cat(2,All(ind).out.exp.rois{:}));
    catch
        rtgs = unique(cat(1,All(ind).out.exp.rois{:}));
    end
    htgs = All(ind).out.exp.targetedCells(rtgs);
    htgs(isnan(htgs))=[];
    cellsToUse = htgs;
    
    %data to build PCA
    %build new composite zdfData;
    dat1 = All(ind).out.vis.zdfData;
    dat2 = All(ind).out.exp.zdfData;
    
    sz1 = size(dat1);
    sz2 = size(dat2);
    mnFrames = min(sz1(2),sz2(2));
    
    
    fullDat = cat(3,All(ind).out.vis.allData(:,1:mnFrames,:),All(ind).out.exp.allData(:,1:mnFrames,:));
    [fulldfDat, fullzdfDat] = computeDFFwithMovingBaseline(fullDat);
    
    dat1 = fullzdfDat(:,:,1:sz1(3));
    dat2 = fullzdfDat(:,:,sz1(3)+1:end);
    
    
    dat1 = dat1(:,1:18,:);
    dat2 = dat2(:,1:18,:);
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
    
    datToUse=datToUse(:,6:18,:);
    
    sz = size(datToUse);
    datToUseMakePCA = reshape(datToUse,[sz(1) sz(2)*sz(3)]);
    datToUseMakePCA = smoothdata(datToUseMakePCA,2); %smooth
    
    [coeff,score,latent,tsquared,explained,mu] = pca(datToUseMakePCA');
    
    
    %data to make vis Epoch
    datToUse = dat1; %All(ind).out.vis.zdfData; %All(ind).out.vis.allData;
    
    datToUse = datToUse(cellsToUse,:,:); %All(ind).out.vis.allData;
    bData = mean(datToUse(:,1:5,:),2);
    datToUse = datToUse-bData; %baseline
    datToUse = smoothdata(datToUse,2); %smooth
    
    sz = size(datToUse);
    datToUse2 = reshape(datToUse,[sz(1) sz(2)*sz(3)]);
    
    % score of vis Data
    visScore = coeff' * datToUse2;
    visScore = reshape(visScore,[sz(1) sz(2) sz(3)]);
    
    colors = colorMapPicker(4,'plasma');
    colors = cat(2,{rgb('grey')},colors,colors);
    
    
    %data to make holo epoch
    hdatToUse = dat2; % All(ind).out.exp.zdfData;
    
    hdatToUse = hdatToUse(cellsToUse,:,:); %All(ind).out.vis.allData;
    bhData = mean(hdatToUse(:,1:5,:),2);
    hdatToUse = hdatToUse-bhData; %baseline
    hdatToUse = smoothdata(hdatToUse,2); %smooth
    
    sz = size(hdatToUse);
    hdatToUse2 = reshape(hdatToUse,[sz(1) sz(2)*sz(3)]);
    
    newHScore = coeff' * hdatToUse2;
    % newScore = reshape(score,[sz(2) sz(3) sz(1)]);
    newHScore = reshape(newHScore,[sz(1) sz(2) sz(3)]);
    
    %define Axis to plot
    ax1 = 1;
    ax2 = 2;
    ax3 = 3;
    
    plotMerge = 1;
    
    figure(3);clf
    hold on
    datToPlot = squeeze(mean(visScore(:,:,visID==uv(1) & trialsToUseVis),3));
    p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
    p.Color = colors{1};
    p.LineWidth = 2;
    
    %plot vis
    for i=1:9;%[1 3 5 7 9];%[1 8 6 4 2];%1:9;%numel(uv)
        
        datToPlot = squeeze(mean(visScore(:,:,visID==uv(i)& trialsToUseVis),3));
        
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
                datToPlot = squeeze(mean(visScore(:,:,(visID==uv(i) | visID==uv(i+4)) & trialsToUseVis),3));
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
    us = All(ind).out.exp.uniqueStims;
    
    stimID = All(ind).out.exp.stimID;
    [s sidx] = sort(All(ind).out.exp.stimParams.Seq);
    us = us(sidx);
    
    if numel(us)<9
        us = [us 0 0 0 0 0 0 0 0 0]; %trailing zeros shouldn't matter
    end
    

    holoIDS = [1 1:9];
    holoPlot = [1 0 1 1 1 1 1 1 1 1 ];

    holoIDS = hIDS{ind};
    holoPlot = hPlot{ind};
    
    plotHoloMerge =1;
    
%     datToPlot = squeeze(mean(newHScore(:,:,stimID==us(2) & trialsToUseExp),3));
%     p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
%     p.Color = colors{1}; %rgb('red'); %colors{i};
%     p.LineWidth =2;
%     p.LineStyle = ':';
    
    
for i=1:numel(us)
    if us(i)~=0
        if ~plotHoloMerge
            datToPlot = squeeze(mean(newHScore(:,:,stimID==us(i) & trialsToUseExp),3));
            p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
            p.Color = colors{holoIDS(i)}; %rgb('red'); %colors{i};
            p.LineWidth =2;
            p.LineStyle = ':';
        else
            if ~holoPlot(i)
                %do nothing
            elseif holoIDS(i)>1 & holoIDS(i)<6; %i>2 && i<7
                newIDS = holoIDS(i)+4;
                tempID = us(find(holoIDS==newIDS));
                if isempty(tempID)
                    tempID=NaN;
                end
                datToPlot = squeeze(mean(newHScore(:,:,(stimID==us(i) | stimID==tempID)  & trialsToUseExp),3));
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
                datToPlot = squeeze(mean(newHScore(:,:,(stimID==us(i) | stimID==tempID)  & trialsToUseExp),3));
                p = plot3(datToPlot(ax1,:),datToPlot(ax2,:),datToPlot(ax3,:));
                p.Color = colors{holoIDS(i)}; %rgb('red'); %colors{i};
                p.LineWidth =2;
                p.LineStyle = ':';
            elseif holoIDS(i)==1
                datToPlot = squeeze(mean(newHScore(:,:,stimID==us(i) & trialsToUseExp),3));
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
    
    bData = mean(dataToUse(:,1:5,:),2);
    dataToUse = dataToUse-bData; %baseline
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
    size(cosSimOut)
    
    testSim = cosSim(2:9,2:9);
    meanSimOnDiag = nanmean(testSim(logical(eye(8))));
    meanSimOffDiag = nanmean(testSim(~logical(eye(8))));
    meanSimOnDiag/meanSimOffDiag;
    
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
     
     for i=2:numel(us);
         r = find(All(ind).out.mani.stimParams.Seq == i-1) ;
         
                  roiPointer = All(ind).out.mani.stimParams.roi{r};
                  roiPointer = unique([All(ind).out.mani.holoRequest.bigListofSequences{i-1}]);
         if roiPointer==0
             cellsHitPerPattern(i)=0;
         else
             cellsTargeted = unique([All(ind).out.mani.rois{roiPointer}]);           
             cellsHitPerPattern(i) = numel(cellsTargeted);
             
             stimPattern = ([All(ind).out.mani.holoRequest.bigListofSequences{i-1}]);
             totalAPsPerPattern(i) = sum(cellsPerROI(stimPattern));
         end
     end
     
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
         e2 = errorbar(1:9,vMean,vSEM);
         xticks(1:9)
         xticklabels({'gray' 0:45:315})
         
         hold on;
         e2 = errorbar(1:9,hMean,hSEM);
         %          xticks(1:9)
         %          xticklabels({'gray' 0:45:315})
         ylabel('zdF/F')
         xlabel('Stimulus')
         pause
         
     end
     end
end

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