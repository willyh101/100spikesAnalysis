%%
% Recreates the key plots of the paper (averaged over ensembles), and 
% reproducing them using a cell-by-cell analysis 
%%
clear; close all; clc;


%%
% Adds all subfolders to the path
restoredefaultpath;
folder = fileparts(which('grayScreenAnalysis.m')); 
addpath(genpath(folder));
rmpath(folder);

%% Specify data location and loadList

% loadPath = '/Users/gregoryhandy/Research_Local/outputdata1/'; % where ever your files are
loadPath = '/Users/greghandy/Research_Local_v2/'; % where ever your files are

% Options: short (for debugging), all (all avaliable outfiles), used (data
% files currently being used here; created for speed)
loadList_all = oriLoadList_GH('used');

%% Loop over each experiment
for exp_loop = 1:length(loadList_all) %10 for 60 offTarget
    
    %% Load the data from a specific experiment (code from Will)
    loadList = loadList_all(exp_loop);
    [All,opts,outVars] = loadData_GH(loadList,loadPath);
           
    %% Additional pre-analysis: 
    % Remove cells from trials where they were offTargets/targets within some
    % buffer zone
    % Key parameters: length of buffer zone and whether to remove targets
    % or offTargets
    % If the total trial count for these cells drops below 10, they are
    % removed entirely
    
    % length(All.out.anal.targets) % total unique targets
    unID = unique(All.out.exp.stimID);
    % Note: The lowest unID is always a gray screen
    % Shift necessary to map onto holoIndex
    IDshift = unID(2)-1;
    totalCells = size(All.out.exp.dfData,1);
    totalTrials = size(All.out.exp.stimID,2);
    
    % Everyone starts being avaliable for analysis
    bufferZone = 1;
    lastTimeStimmed = bufferZone+zeros(totalCells,1);
    % Loop through all trials sequentially
    for ii = 1:totalTrials
        % Remove these cells from this trial
        doNotInclude = lastTimeStimmed<bufferZone;
        All.out.exp.dfData(doNotInclude,:,ii) = nan;
        
        % Update the lastTimeStimmed list
        % Figure out who is currently stimmed
        holoIndex = All.out.exp.stimID(ii)-IDshift;
        if holoIndex > 0
            % If you only want to remove targets, uncomment this line
            %currTargets = All.out.exp.holoTargets{holoIndex}(~isnan(All.out.exp.holoTargets{holoIndex}));
            currTargets=find(All.out.anal.offTargetRisk(holoIndex,:)==1);
        else
            % no stim, so no current targets 
            currTargets=[];
        end
        lastTimeStimmed = lastTimeStimmed + 1;
        % Reset the count of those just stimmed
        lastTimeStimmed(currTargets) = 0;
    end
    
    %% Cells to use
    ROIinArtifact  = All.out.anal.ROIinArtifact; % cells at edge of FOV in laser artifact
    cellsToUse = ~ROIinArtifact' & ~All.out.anal.cellsToExclude;
       
    %% Control response (using the no stim epochs)
    % Method 1:
    stimIDs = unique(All.out.exp.stimID);
    trialsToUse = ismember(All.out.exp.stimID, stimIDs(1)) &...
        All.out.exp.lowMotionTrials &...
        All.out.exp.lowRunTrials & All.out.exp.visID==1;
    tempTotalTrials = sum(trialsToUse);
    
    recWinSec=opts.recWinRange+All.out.exp.visStart;
    winToUse = round(recWinSec*All.out.info.FR);
    bwinToUse = max(floor([0 All.out.exp.visStart]*All.out.info.FR),[1 1]);
    tracesHolodfData = All.out.exp.dfData(cellsToUse, :, trialsToUse);
    dfCellResp = nanmean((nanmean(tracesHolodfData(:,winToUse(1):winToUse(end),:),2)),3);
    baselineEst = nanmean((nanmean(tracesHolodfData(:,bwinToUse(1):bwinToUse(2),:),2)),3);
    dffCellCtlResp_noStim = (dfCellResp-baselineEst);
    
    % Remove cells that don't have enough (10) trials
    badCells = (tempTotalTrials-sum(isnan(tracesHolodfData(:,1,:)),3))<10;
    dffCellCtlResp_noStim(badCells) = nan;
    
    % Method 2: (used as a check for bufferZone=0; doesn't work for bufferZone>0)
    tempCtlResp = squeeze(All.out.anal.respMat(1,cellsToUse));
    tempCtlBase = squeeze(All.out.anal.baseMat(1,cellsToUse));
    dffCellCtlResp_v3 = (tempCtlResp-tempCtlBase)';
    if bufferZone<=0 && norm(dffCellCtlResp_noStim-dffCellCtlResp_v3)~=0
       error('Something went wrong calculating the gray screen response')
    end
    
    figure(100+exp_loop); clf;
    subplot(1,2,1); hold on;
    tempAve = nanmean(tracesHolodfData,3)-baselineEst;
    tempAve(badCells,:) = [];
    tracesHolodfDataTemp = tracesHolodfData(~badCells,:,:);
    bdataTemp = All.out.exp.bdata(~badCells,trialsToUse);
    bdataTemp2 =  All.out.anal.baseMat(:,~badCells);
    [vals,Indices]=sort(dffCellCtlResp_noStim(~badCells),'ascend');
    
    imagesc(tempAve(Indices,:))
    colorbar
    colormap rdbu
    caxis([-0.5 0.5])
    set(gca,'fontsize',16)
    xticks([1:1:5]*All.out.info.FR)
    xticklabels([1 2 3 4])
    ylim([0 length(Indices)])
    xlim([1 size(tempAve,2)])
    title(loadList{1},'Interpreter', 'none')
    xlabel('Time (sec)')
    
    timeVecTemp = [0:1:size(tempAve,2)-1]/All.out.info.FR;    
    subplot(1,2,2); hold on;
    tempMean = nanmean(tempAve);
    tempSTDErr = nanstd(tempAve)/sqrt(size(tempAve,1));
    plot(timeVecTemp,mean(tempAve),'linewidth',1.5)
    hTemp = patch([timeVecTemp fliplr(timeVecTemp)], [tempMean+tempSTDErr fliplr(tempMean-tempSTDErr)],[0 0.447 0.741]);
    hTemp.FaceAlpha=0.2;
    hTemp.EdgeColor=[1 1 1];
    xlim([0 timeVecTemp(end)])
    ylim([-0.1 0.1])
    ylimTemp = ylim;
    xlimTemp = xlim;
    plot(xlimTemp,0*xlimTemp,'k--')
    plot(All.out.exp.visStart+0*ylimTemp,ylimTemp,'k--')
    plot(All.out.exp.visStart+1+0*ylimTemp,ylimTemp,'k--')
    set(gca,'fontsize',16)
    title(sprintf('Population resp: %.3f',nanmean(vals)))
    xlabel('Time (sec)')
    
    if unique(All.out.anal.minDistbyHolo(1,:)) ~=0
        error('Something off with minDistbyHolo');
    end
        
    fprintf('Exp %d out of %d done loading\n',exp_loop, length(loadList_all));
end

