function [All outVars] = CalcPVisRFromGMN(All,opts,outVars)
%%
try 
    recWinRange = opts.visRecWinRange;
catch
    disp('Might not have the visRecWinRange using recWinRange Instead')
    recWinRange =opts.recWinRange;%
end
visAlpha = opts.visAlpha;
numExp = numel(All);

for ind = 1:numExp;
    if ~isfield(All(ind).out,'vis2')
        disp(['No Vis2 Field ind: ' num2str(ind)]);
        visPercent(ind) = 0; 
    else
        visStart = All(ind).out.vis2.visStart;
        recWinSec = recWinRange + visStart; %recording window relative to when vis start
        sz = size(All(ind).out.vis2.zdfData);
        
        winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
        bwinToUse = max(round([0 All(ind).out.vis2.visStart]*All(ind).out.info.FR),[1 1]);
        All(ind).out.vis2.win = winToUse;
        All(ind).out.vis2.bwin = bwinToUse;
        
        vrdata = squeeze(mean(All(ind).out.vis2.zdfData(:,winToUse(1):winToUse(2),:),2));
        vbdata = squeeze(mean(All(ind).out.vis2.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
        All(ind).out.vis2.rdData=vrdata;
        All(ind).out.vis2.bdata=vbdata;
        
        %ID tuned Cells, should comparing no contrast to with contrast
        pVisR=[];%pVisT=[];
        for i=1:All(ind).out.anal.numCells
            trialsToUse = All(ind).out.vis2.visID~=0;% &...
%                 All(ind).out.vis2.lowMotionTrials &...
%                 All(ind).out.vis2.lowRunTrials;% &...
            %             All(ind).out.vis2.stimID==min(All(ind).out.vis2.stimID);
            if size(trialsToUse,2) > size(All(ind).out.vis2.rdata,2)
                trialsToUse = trialsToUse(1:size(All(ind).out.vis2.rdata,2));
            end
            
            
            pVisR(i) = anova1(All(ind).out.vis2.rdData(i,trialsToUse),All(ind).out.vis2.visID(trialsToUse),'off');
            %         pVisT(i) = anova1(All(ind).out.vis2.rdata(i,trialsToUse),All(ind).out.vis2.visID(trialsToUse),'off');
        end
        All(ind).out.anal.pVisR = pVisR;
        All(ind).out.anal.visPercent = sum(pVisR<visAlpha) / numel(pVisR);
        visPercent(ind) =  All(ind).out.anal.visPercent;
        
        outVars.pVisR{ind} = pVisR;
    end
end
outVars.visPercent = visPercent;
disp('Calculated pVisR from out.vis')
