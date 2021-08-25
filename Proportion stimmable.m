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
        if numel(All(ind).out.stm.stimID) ~= size(All(ind).out.stm.zdfData,3)
            disp(['Trimming Ind: ' num2str(ind)]);
            All(ind).out.stm.stimID = All(ind).out.stm.stimID(3:end);
            All(ind).out.stm.runVal = All(ind).out.stm.runVal(3:end,:);
        end
        
    elseif strcmp([All(ind).out.info.date '_' All(ind).out.info.mouse],'200810_i140_2')
        if numel(All(ind).out.stm.stimID) ~= size(All(ind).out.stm.zdfData,3)
            disp(['Trimming Ind: ' num2str(ind)]);
            All(ind).out.stm.stimID = All(ind).out.stm.stimID(1:61);
            All(ind).out.stm.runVal = All(ind).out.stm.runVal(1:61,:);
        end
        
    elseif strcmp([All(ind).out.info.date '_' All(ind).out.info.mouse],'200901_w19_1')
       
            disp(['Correcting Ind: ' num2str(ind)]);
            All(ind).out.stm.stimID(All(ind).out.stm.stimID==7)=0;
            
        
    end
   
end

%% Identify the Experiment type for comparison or exclusion
outVars=[];
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;

%% Exclude List

exclList = [2 6 8 10 17];[6 17 8 10];%[6 8 10] ;

%% run
        figure(8);clf; hold on
        figure(7);clf
        
        summaryDataPowers = 0:25:150;
        summaryData = nan([numel(summaryDataPowers) numExps]);
        isAi203=nan([1 numExps]);
for ind = 1:numExps
    if isfield(All(ind).out, 'stm') && ~ismember(ind,exclList)
        stm = All(ind).out.stm;
        info = All(ind).out.info;
        
        if ~isfield(info,'FR')
            info.FR = 6;
        end
        
        dataToUse = stm.zdfData;
        runVector =stm.runVal;
        targettedCells =stm.targetedCells;
        try
        outputPowers = cellfun(@(x) str2num(x(regexp(x,'\. ','once')+2:regexp(x,'mW','once')-1)),stm.outputsInfo.OutputNames,'uniformoutput',1);
        catch
            disp(['name error. ind ' num2str(ind)])
            outputPowers=[];
        end
            
          try
              outputStimID = stm.outputsInfo.OutputStimID;
          catch
            outputStimID= stm.uniqueStims;
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
            
            if ~isempty(timeOfStim) & timeOfStim*info.FR<size(dataToUse,2)-0.5*info.FR
                frameOfStim = round(timeOfStim*info.FR);
                winRS = frameOfStim:min(frameOfStim+3-1,size(runVector,2));
                runSpeed = runVector(:,winRS)';
                
                runTrialsToInclude = mean(runSpeed)<10000.5;
                if numel(runTrialsToInclude)>size(dataToUse,3)
                    runTrialsToInclude = runTrialsToInclude(1:size(dataToUse,3));
                end
                
                temp =  dataToUse(i,frameOfStim:min(frameOfStim+FramesToCountStim-1,size(dataToUse,2)),runTrialsToInclude);
                if size(temp,2)<FramesToCountStim
                    padSize = FramesToCountStim - size(temp,2);
                    temp = padarray(temp,[0 padSize 0],nan,'post');
                end
                stimTestResp(i,:,runTrialsToInclude) = temp;
                
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
        
        
        %%Proportion stimmable
        stimID = stm.stimID;
%         if numel(stimID)>size(dataToUse,3)
%             stimID = stimID(3:end);%size(dataToUse,3));
%         end
        prob1=[];
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
        
        fProb1 = double(prob1<0.05);
        fProb1(isnan(prob1))=NaN;
        
        fProb1Sort = fProb1(:,sIdx);
        figure(8);
        if outVars.indExpressionType(ind)==5
        plot(powersSorted,nanmean(fProb1Sort),'color',[rgb('DarkOrange') 0.25],'lineWidth',1)
        isAi203(ind)=1;
        else
            plot(powersSorted,nanmean(fProb1Sort),'color',[rgb('ForestGreen') 0.25],'lineWidth',1)
                    isAi203(ind)=0;
        end

        
        figure(7)
        subplot(6,5,ind)
        if outVars.indExpressionType(ind)==5
            plot(powersSorted,mean(fProb1Sort),'color',[rgb('DarkOrange') 1],'lineWidth',1);
        else
            plot(powersSorted,mean(fProb1Sort),'color',[rgb('ForestGreen') 1],'lineWidth',1)
        end
        xlim([0 175]);
        ylim([0 1]);
        title({['ind: ' num2str(ind)]; [info.date '_' info.mouse]; [num2str(numel(testedCells)) ' Cells']},'Interpreter','none')
        drawnow
        %         pause
        temp=nan([1 numel(summaryDataPowers)]);
        avData = mean(fProb1Sort);
        for k=1:numel(powersSorted)
            f = find(powersSorted(k)==summaryDataPowers);
            temp(f)= avData(k); 
            if isempty(f)
                disp(['No Match ind: ' num2str(ind) ' pwr: ' num2str(powersSorted(k))])
            end
        end
        summaryData(:,ind)=temp;
%         
    else
        disp(['Did not run ind: ' num2str(ind) ' ' num2str(All(ind).out.info.date) '_' All(ind).out.info.mouse])
    end
end

        figure(8);

        
        xlim([0 155]);
        ylim([0 1]);
        xlabel('Power (mW)')
        ylabel('Proportion Stimmable')
        xticks(0:25:150)
        
        disp('done')
        sumDat1 =  summaryData(:,isAi203==1);
        mn1 = nanmean(sumDat1,2);
        std1 = nanstd(sumDat1,[],2);
        npts1 = sum(~isnan(sumDat1)');
        ste1 = std1./sqrt(npts1');
        
        p1 = errorbar(summaryDataPowers,mn1,ste1);
        
%         p = plot(summaryDataPowers, nanmean(summaryData(:,isAi203==1),2));
        p1.Color = rgb('DarkOrange');
        p1.LineWidth = 2;
        p1.Marker='o';
        p1.MarkerFaceColor = rgb('DarkOrange');
        p1.CapSize=0
        
        sumDat2 = summaryData(:,isAi203==0);
        mn2 = nanmean(sumDat2,2);
        std2 = nanstd(sumDat2,[],2);
        npts2 = sum(~isnan(sumDat2)');
        ste2 = std2./sqrt(npts2');
        
        p1 = errorbar(summaryDataPowers,mn2,ste2);
        
        p1.Color = rgb('ForestGreen');
        p1.LineWidth = 2;
        p1.Marker='o';
        p1.MarkerFaceColor = rgb('ForestGreen');
        p1.CapSize=0;
        
        set(gca,'TickDir','out')
%% Repeated Measures Anova
sData = summaryData';

sData(isnan(isAi203),:)= [];
isAi = isAi203;
isAi(isnan(isAi203))=[];

t = table(isAi',sData(:,1),sData(:,2),sData(:,3),sData(:,4),sData(:,5),sData(:,6),sData(:,7),...
    'VariableNames',{'isAi203','p0','p25','p50','p75','p100','p125','p150'});


    rm = fitrm(t,'p0-p125 ~ isAi203','WithinDesign',summaryDataPowers(1:6));
    
    ranova(rm)