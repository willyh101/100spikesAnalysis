numExps = numel(All);
ensemblesToUse = outVars.ensemblesToUse; 
IndsUsed = outVars.IndsUsed;
minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));

ind=2;

us = unique(All(ind).out.exp.stimID);
vs = unique(All(ind).out.exp.visID);
vs(vs==0)=[];

ROIinArtifact = All(ind).out.anal.ROIinArtifact;
offTargetRisk = All(ind).out.anal.offTargetRisk;

strtFrame = All(ind).out.anal.recStartFrame;
newStart = strtFrame-minStrtFrame+1;

pVisR = All(ind).out.anal.pVisR;
osis = All(ind).out.anal.osi;

trialsToUse = All(ind).out.anal.defaultTrialsToUse;

i = randperm(numel(us),1);
s = us(i);

h = All(ind).out.exp.stimParams.roi{i};

if h>0
    tg = All(ind).out.exp.holoTargets{h};
    tg(isnan(tg))=[];
else
    tg=[];
end
cellList = 1:numel(ROIinArtifact);

if i==1
    cellsToUse = ~ROIinArtifact' & pVisR<0.05 ;
    cellIdxs = find(cellsToUse);
else
    cellsToUse = ~ROIinArtifact' &...
        ~offTargetRisk(h,:) &...
        ~ismember(cellList,tg) &...
        pVisR<0.05 &...
        osis>0.5;
    cellIdxs = find(cellsToUse);
end

dat = All(ind).out.exp.zdfData(cellsToUse,newStart:end,trialsToUse &...
    All(ind).out.exp.stimID==s &...
    1);

mDat = mean(dat,3);
mmDat = mean(mDat,1); %pop Average
sdDat = std(mDat);
nDat = size(mDat,1);


figure(666);
clf
% size(mDat)
cell2plot = randperm(size(mDat,1),1);
data = mDat(cell2plot,:);
% err = 
plot(data);
            

