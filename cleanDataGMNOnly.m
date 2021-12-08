function [All] = cleanDataGMNOnly(All,opts)
%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=opts.FRDefault;%6;
recWinRange =opts.recWinRange;% [0.5 1.5];% %from vis Start [1.25 2.5];


%Stim Success Thresholds
stimsuccessZ = opts.stimsuccessZ;%0.25; %over this number is a succesfull stim
stimEnsSuccess = opts.stimEnsSuccess;% 0.5; %fraction of ensemble that needs to be succsfull

%run Threshold
runThreshold = opts.runThreshold;% 6 ; %trials with runspeed below this will be excluded

numExps = numel(All);

clear ensStimScore

for ind =1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    if ~isfield(All(ind).out.info,'FR')
        All(ind).out.info.FR=FRDefault;
    end
    All(ind).out.anal.numCells = size(All(ind).out.exp.zdfData,1);
    %Vis Section
    if isfield(All(ind).out,'vis2')
        sz2 = size(All(ind).out.vis2.zdfData);
        try
            visStart = All(ind).out.vis2.visStart;
        catch
            visStart = 0.92; %seemed typical
            All(ind).out.vis2.visStart = visStart;
            fprintf('\nError vis2.visStart not detected...')
        end
        winToUse = min(round((visStart+recWinRange).*All(ind).out.info.FR),[inf sz2(2)]) ;
        bwinToUse = max(round([0 visStart]*All(ind).out.info.FR),[1 1]);
        %      winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz2(2)]) ;
        %      bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);
        All(ind).out.anal.visRecWinUsed = winToUse;
        
        
        rdata = squeeze(mean(All(ind).out.vis2.zdfData(:,winToUse(1):winToUse(2),:),2));
        bdata = squeeze(mean(All(ind).out.vis2.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
        
        All(ind).out.vis2.rdata=rdata;
        All(ind).out.vis2.bdata=bdata;
        
        if ~isfield(All(ind).out.vis,'lowRunTrials')
            try
                runVal = All(ind).out.vis2.runVal;
                rnSz = size(runVal);
                runperiod = [1:min(All(ind).out.anal.visRecWinUsed(2),rnSz(2))];
                lowRunVals = mean((runVal(:,runperiod)<runThreshold)');
                lowRunTrials = lowRunVals>0.75; %percent of frames that need to be below run threshold
                All(ind).out.vis2.lowRunTrials = lowRunTrials;
            catch
                disp('Vis Run Trial Problem')
            end
        end
        fprintf('\n')
    end
end
