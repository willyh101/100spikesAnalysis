function dat = makeMeanResponse(out, stimID, visID, baseline)

if isfield(out.exp, 'dataToUse')
    dataToUse = out.exp.dataToUse;
else
    disp(['ind ' num2str(ind) '. no data to use, using zdfData']);
    dataToUse = out.exp.zdfData;
end

ROIinArtifact = out.anal.ROIinArtifact;
offTargetRisk = out.anal.offTargetRisk;

newStart = out.anal.recStartFrame;

pVisR = out.anal.pVisR;

trialsToUse = out.anal.defaultTrialsToUse;

us = unique(out.exp.stimID);
vs = unique(out.exp.visID);

if ~ismember(stimID, us)
    error(['StimID of ' num2str(stimID) ' not in range of unique stims.'])
end
if ~ismember(visID, vs)
    error(['visID of ' num2str(visID) ' not in range of unique vis stims.'])
end

s = stimID;
v = visID;
i = find(us == s);
h = out.exp.stimParams.roi{i};

if h>0
    tg = out.exp.holoTargets{h};
    tg(isnan(tg))=[];
else
    tg=[];
end

cellList = 1:numel(ROIinArtifact);

if h==0 %changed to account for when things don't run in order
    cellsToUse = ~ROIinArtifact' & out.anal.cellsToInclude;% & pVisR<0.05 ;
else
    cellsToUse = ~ROIinArtifact' &...
        ~offTargetRisk(h,:) &...
        ~ismember(cellList,tg) &...
        out.anal.cellsToInclude;% & pVisR<0.05;
end

dat = dataToUse(cellsToUse,:,trialsToUse &...
                out.exp.stimID==s &...
                out.exp.visID==v );

if baseline
    dat = dat-mean(dat(:,1:newStart),2);
end

dat = squeeze(mean(dat, 3));