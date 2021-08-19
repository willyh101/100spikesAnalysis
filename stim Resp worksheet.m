stm;
info;




dataToUse = stm.zdfData;
runVector =stm.runVal;

outputPowers = cellfun(@(x) str2num(x(regexp(x,'\. ','once')+2:regexp(x,'mW','once')-1)),stm.outputsInfo.OutputNames,'uniformoutput',1);
outputStimID = stm.outputsInfo.OutputStimID;

fst = stm.holoRequest.bigListOfFirstStimTimes(:,1); %first stim times

numCells = size(dataToUse,1);
numTrials = size(dataToUse,3);
FramesToCountStim = 8; 
FramesToCountBase =4;


stimTestResp = nan([numCells, FramesToCountStim,size(dataToUse,3)]);
preTestResp = nan([numCells, FramesToCount,size(dataToUse,3)]);
midVal = nan([numCells, size(dataToUse,3)]);
testedCells =[];

for i=1:numCells
    tgNum = find(targettedCells==i,1);
    timeOfStim = fst(tgNum);
    timeOfStim(isnan(timeOfStim)) =[];
    
    if ~isempty(timeOfStim)
        frameOfStim = round(timeOfStim*info.FR);
        winRS = frameOfStim:min(frameOfStim+3-1,size(runVector,2));
        runSpeed = runVector(:,winRS)';
        
        runTrialsToInclude = mean(runSpeed)<0.5;
        
        stimTestResp(i,:,runTrialsToInclude) = dataToUse(i,frameOfStim:frameOfStim+FramesToCountStim-1,runTrialsToInclude);
        preTestResp(i,:,runTrialsToInclude) = dataToUse(i,max(1,frameOfStim-FramesToCountBase-1):frameOfStim-2,runTrialsToInclude);
        midVal(i,runTrialsToInclude) = dataToUse(i,frameOfStim-1,runTrialsToInclude); 

        testedCells(end+1)=i;
    end
    
    
end

[powersSorted sIdx] = sort(outputPowers);

%% heatmap
meanTraces=[];
for i=1:numel(testedCells)
    c = testedCells(i);
    for k = 1:numel(outputStimID)
        s = outputStimID(k);        
        datMat = cat(1,squeeze(preTestResp(c,:,stimID==s)),midVal(c,stimID==s),squeeze(stimTestResp(c,:,stimID==s)));
        meanTrace = nanmean(datMat');
        
%         baseVal = nanmean(nanmean(preTestResp(c,:,stimID==s)));
%         meanTrace = meanTrace - baseVal;
        
        meanTraces(i,:,k) = meanTrace;
        
%          datToPlot(i,k,:) = 
    end
end
figure(5);clf
for i =1:numel(outputStimID)
    p = sIdx(i);
    subplot(3,5,i)
    imagesc(meanTraces(:,:,p))
    title(outputPowers(p))
    caxis([0 3])
end

%% Proportion stimmable
for i =1:numel(testedCells)
    c = testedCells(i);
    fprintf('. ');
    for k = 1:numel(outputStimID)
        s = outputStimID(k);
        baseVals = mean(squeeze(preTestResp(c,:,stimID==s)));
        stimVals = mean(squeeze(stimTestResp(c,:,stimID==s)));
        
        prob1(i,k) = ranksum(stimVals,baseVals,'tail','right');
    end
end

fProb1 = prob1<0.05;

fProb1Sort = fProb1(:,sIdx);
figure(6);clf
plot(powersSorted,mean(fProb1Sort),'color','k','lineWidth',2)
%% ranksum to compare
prob = nan([numel(testedCells) numTrials]);
for i =1:numel(testedCells)
    c = testedCells(i);
    fprintf('. ');
    for k = 1:numTrials
        if all(isnan([stimTestResp(c,:,k) preTestResp(c,:,k)]))
            prob(i,k) = nan;
        else
            prob(i,k) = ranksum(stimTestResp(c,:,k),preTestResp(c,:,k),'tail','right');
        end
    end
end
disp done;

fProb=[];
for i =1:numel(outputStimID)    
    S =outputStimID(i);
    
    fProb(:,i) = nanmean(prob(:,stimID==S)<0.05,2);
end

    
fProbSort = fProb(:,sIdx);
figure(2);clf
plot(powersSorted,fProbSort')
hold on; plot(powersSorted,mean(fProbSort),'color','k','lineWidth',2)

%% based on diff from 0 cond
hitCount = nan([numel(testedCells) numTrials]);

s0 =outputStimID(find(outputPowers==0));
zeroData = dataToUse(:,:,stimID==s0);
sz = size(zeroData);
zeroData = reshape(zeroData,[sz(1) sz(2)*sz(3)]);

m0Data = nanmean(zeroData');
s0Data = nanstd(zeroData');

stdThresh = 2;


for i =1:numel(testedCells)
    c = testedCells(i);
    
    for k = 1:numTrials
        hitCount(i,k) = any ( stimTestResp(c,:,k) > m0Data(c)+stdThresh*s0Data(c));
    end
end

for i =1:numel(outputStimID)    
    S =outputStimID(i); 
    fProb(:,i) = mean(hitCount(:,stimID==S),2);
end
[s sIdx] = sort(outputPowers);
    
fProbSort = fProb(:,sIdx);
figure(3);clf
plot(powersSorted,fProbSort')
hold on; plot(powersSorted,mean(fProbSort),'color','k','lineWidth',2)


%% Binning in 1s increments on diff from 0 cond
hitCount = nan([numel(testedCells) numTrials]);

s0 =outputStimID(find(outputPowers==0));
zeroData = dataToUse(:,:,stimID==s0);
sz = size(zeroData);
zeroData = reshape(zeroData,[sz(1) sz(2)*sz(3)]);

sz= size(zeroData);
clear zeroData2
for i = 1:floor(sz(2)/FramesToCountStim)
    zeroData2(:,i) = nanmean(zeroData(:,(i-1)*FramesToCountStim+1:i*FramesToCountStim),2);
end

m0Data = nanmean(zeroData2');
s0Data = nanstd(zeroData2');

stdThresh = 2;


for i =1:numel(testedCells)
    c = testedCells(i);
    
    for k = 1:numTrials
        hitCount(i,k) =  mean(stimTestResp(c,:,k)) > m0Data(c)+stdThresh*s0Data(c);
    end
end

for i =1:numel(outputStimID)    
    S =outputStimID(i); 
    fProb(:,i) = mean(hitCount(:,stimID==S),2);
end
[s sIdx] = sort(outputPowers);
    
fProbSort = fProb(:,sIdx);
figure(4);clf
plot(powersSorted,fProbSort')
hold on; plot(powersSorted,mean(fProbSort),'color','k','lineWidth',2)

%% Proportion stimmable
