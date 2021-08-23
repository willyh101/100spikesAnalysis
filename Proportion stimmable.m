ai203LoadList;
% testLoadList
loadPath = 'T:\Outfiles';

%% load
numExps = numel(loadList);

clear All
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end
disp('Loaded')

%% error check
for ind = 1:numExps
    if strcmp([All(ind).out.info.date '_' All(ind).out.info.mouse],'200727_w26_1')
        All(ind).out.stm.stimID = All(ind).out.stm.stimID(3:end);
        All(ind).out.stm.runVal = All(ind).out.stm.runVal(3:end,:);
    end
end

%% Identify the Experiment type for comparison or exclusion
outVars=[];
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;

%% run

for ind = 1:numExps
    if isfield(All(ind).out, 'stm')
        stm = All(ind).out.stm;
        
        dataToUse = stm.zdfData;
        runVector =stm.runVal;
        targettedCells =stm.targetedCells;
        try
        outputPowers = cellfun(@(x) str2num(x(regexp(x,'\. ','once')+2:regexp(x,'mW','once')-1)),stm.outputsInfo.OutputNames,'uniformoutput',1);
        outputStimID = stm.outputsInfo.OutputStimID;
        catch
            disp('name error')
        end
        
        fst = stm.holoRequest.bigListOfFirstStimTimes(:,1); %first stim times
        
        numCells = size(dataToUse,1);
        numTrials = size(dataToUse,3);
        FramesToCountStim = 8;
        FramesToCountBase =4;
        
        
        stimTestResp = nan([numCells, FramesToCountStim,size(dataToUse,3)]);
        preTestResp = nan([numCells, FramesToCountBase,size(dataToUse,3)]);
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
                if numel(runTrialsToInclude)>size(dataToUse,3)
                    runTrialsToInclude = runTrialsToInclude(1:size(dataToUse,3));
                end
                
                stimTestResp(i,:,runTrialsToInclude) = dataToUse(i,frameOfStim:frameOfStim+FramesToCountStim-1,runTrialsToInclude);
                temp = dataToUse(i,max(1,frameOfStim-FramesToCountBase-1):frameOfStim-2,runTrialsToInclude);
                
                if size(temp,2)<FramesToCountBase
                    padSize = FramesToCountBase - size(temp,2);
                    temp = padarray(temp,[0 padSize 0],nan,'pre');
                end
                
                preTestResp(i,:,runTrialsToInclude) =temp;
                midVal(i,runTrialsToInclude) = dataToUse(i,frameOfStim-1,runTrialsToInclude);
                
                testedCells(end+1)=i;
            end
            
            
        end
        
        [powersSorted sIdx] = sort(outputPowers);
        
        
        %% Proportion stimmable
        stimID = stm.stimID;
%         if numel(stimID)>size(dataToUse,3)
%             stimID = stimID(3:end);%size(dataToUse,3));
%         end
        
        for i =1:numel(testedCells)
            c = testedCells(i);
            fprintf('. ');
            for k = 1:numel(outputStimID)
                s = outputStimID(k);
                baseVals = nanmean(squeeze(preTestResp(c,:,stimID==s)));
                stimVals = nanmean(squeeze(stimTestResp(c,:,stimID==s)));
                
                try
                    prob1(i,k) = ranksum(stimVals,baseVals,'tail','right');
                catch
                    prob1(i,k) = nan;
                end
            end
        end
        disp('.')
        
        fProb1 = prob1<0.05;
        
        fProb1Sort = fProb1(:,sIdx);
        figure(6);clf
        plot(powersSorted,mean(fProb1Sort),'color','k','lineWidth',2)
        title(num2str(ind))
        pause
        
    end
end
