% just plot one experiment the way I want it


minStrtFrame = min(arrayfun(@(x) x.out.anal.recStartFrame,All));

allMeanTS = outVars.allMeanTS;
allStdTS = outVars.allStdTS;
allnumTS = outVars.allnumTS;

baseline=1;

shortestRec = min(cellfun(@(x) size(x,2),allMeanTS));
temp = cellfun(@(x) x(2:end,1:shortestRec),allMeanTS,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
meanTSSquare = cat(1,temp{:});

%no Resp
temp = cellfun(@(x) x(1,1:shortestRec),allMeanTS,'uniformoutput',0);
if baseline
    temp = cellfun(@(x) x-mean(x(:,1:minStrtFrame),2),temp,'uniformoutput',0);
end
meanTSSquareNR = cat(1,temp{:});

temp = cellfun(@(x) x(2:end,1:shortestRec),allStdTS,'uniformoutput',0);
stdTSSquare = cat(1,temp{:});

%% stuff to plot

ind = 42; % 43 and 42 suspect, 36 looks normal

% for ind = IndsUsed

ensIdx = outVars.ensIndNumber == ind;
ensNum = find(ensIdx);

figure(66)
clf
sgtitle(['Ind = ' num2str(ind)])

subplot(1,4,1)
hold on
plot(meanTSSquareNR(ind, :))
plot(mean(meanTSSquare(ensNum, :)))
legend({'Control', 'Holo'})
title('Mean Responses')
ylabel('\DeltaZ-Scored dF/F')

subplot(1,4,2)
ctrlCond = min(All(ind).out.exp.stimID);
resp = makeMeanResponse(All(ind).out, ctrlCond, 1, 1);
imagesc(resp)
title('Cell Response (control)')
caxis([-2 2])
colorbar

subplot(1,4,3)
hold on
for i = ensNum
    plot(meanTSSquare(i,:))
end
plot(mean(meanTSSquare(ensNum,:)), 'LineWidth',3, 'Color','k')
title('Responses to Holos')
ylabel('\DeltaZ-Scored dF/F')

subplot(1,4,4)
imagesc(meanTSSquare(ensNum,:))
colormap rdbu
caxis([-0.05 0.05])
colorbar
title('Ensemble Respones')
% pause
% end