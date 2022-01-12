function [All] = cleanDataVisOnly(All,opts)

numExps = numel(All);

for ind =1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    if isfield(All(ind).out.exp, 'dataToUse')
        dataToUse = All(ind).out.exp.dataToUse;
    else
        disp(['ind ' num2str(ind) '. no data to use, using zdfData']);
        dataToUse = All(ind).out.exp.zdfData;
    end
    if ~isfield(All(ind).out.info,'FR')
        All(ind).out.info.FR=FRDefault;
    end
    %Vis Section
    if isfield(All(ind).out,'vis')
        %         sz2 = size(All(ind).out.vis.zdfData);
        try
            visStart = All(ind).out.vis.visStart;
        catch
            visStart = 0.92; %seemed typical
            All(ind).out.vis.visStart = visStart;
            fprintf('\nError vis.visStart not detected...')
        end
        %         winToUse = min(round((visStart+recWinRange).*All(ind).out.info.FR),[inf sz2(2)]) ;
        %         bwinToUse = max(round([0 visStart]*All(ind).out.info.FR),[1 1]);
        %         %      winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz2(2)]) ;
        %         %      bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);
        %         All(ind).out.anal.visRecWinUsed = winToUse;
        %
        %
        %         rdata = squeeze(mean(All(ind).out.vis.zdfData(:,winToUse(1):winToUse(2),:),2));
        %         bdata = squeeze(mean(All(ind).out.vis.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
        %
        %         All(ind).out.vis.rdata=rdata;
        %         All(ind).out.vis.bdata=bdata;
        
            sz = size(dataToUse);
        recWinSec = opts.visRecWinRange + All(ind).out.vis.visStart; %recording window relative to when vis start
        
        winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
        bwinToUse = max(floor([0 All(ind).out.vis.visStart]*All(ind).out.info.FR),[1 1]);
        
        All(ind).out.anal.recWinUsed = winToUse;
        All(ind).out.anal.bwinToUse = bwinToUse;
        All(ind).out.anal.recStartTime = visStart;
        All(ind).out.anal.recStartFrame = round(visStart*All(ind).out.info.FR);

        runThreshold = opts.runThreshold;% 6 ; %trials with runspeed below this will be excluded
runValPercent = opts.runValPercent;
        
        if ~isfield(All(ind).out.vis,'lowRunTrials')
            try
                runVal = All(ind).out.vis.runVal;
                rnSz = size(runVal);
                runperiod = [1:min(All(ind).out.anal.recWinUsed(2),rnSz(2))];
                lowRunVals = mean((runVal(:,runperiod)<runThreshold)');
                lowRunTrials = lowRunVals>runValPercent; %percent of frames that need to be below run threshold
                All(ind).out.vis.lowRunTrials = lowRunTrials;
            catch
                disp('Vis Run Trial Problem')
            end
        end
    end
end