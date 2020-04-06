function [All outVars] = CalcPVisRFromVis(All,opts,outVars)
%%
recWinRange =opts.recWinRange;%
visAlpha = opts.visAlpha;
numExp = numel(All);

for ind = 1:numExp;
    if ~isfield(All(ind).out,'vis')
        disp(['No Vis Field ind: ' num2str(ind)]);
        visPercent(ind) = outVars.visPercent(ind); 
    else
        visStart = All(ind).out.vis.visStart;
        recWinSec = recWinRange + visStart; %recording window relative to when vis start
        sz = size(All(ind).out.vis.zdfData);
        
        winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
        bwinToUse = max(round([0 All(ind).out.exp.visStart]*All(ind).out.info.FR),[1 1]);
        
        vrdata = squeeze(mean(All(ind).out.vis.zdfData(:,winToUse(1):winToUse(2),:),2));
        vbdata = squeeze(mean(All(ind).out.vis.zdfData(:,bwinToUse(1):bwinToUse(2),:),2));
        All(ind).out.vis.rdData=vrdata;
        All(ind).out.vis.bdata=vbdata;
        
        %ID tuned Cells, should comparing no contrast to with contrast
        pVisR=[];%pVisT=[];
        for i=1:All(ind).out.anal.numCells
            trialsToUse = All(ind).out.vis.visID~=0 &...
                All(ind).out.vis.lowMotionTrials &...
                All(ind).out.vis.lowRunTrials;% &...
            %             All(ind).out.vis.stimID==min(All(ind).out.vis.stimID);
            
            
            
            pVisR(i) = anova1(All(ind).out.vis.rdData(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
            %         pVisT(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
        end
        All(ind).out.anal.pVisR = pVisR;
        All(ind).out.anal.visPercent = sum(pVisR<visAlpha) / numel(pVisR);
        visPercent(ind) =  All(ind).out.anal.visPercent;
    end
    
end
outVars.visPercent = visPercent;
disp('Calculated pVisR from out.vis')
