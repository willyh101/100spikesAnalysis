%%
% Plots results pertaining to the targeted cells
%%
function FigTargets(cellTable)

totalNumEns = cellTable.ensNum(end);

tgRespAve = zeros(totalNumEns,1); tgRespErr = zeros(totalNumEns,1); 
numTg = zeros(totalNumEns,1); numTgAct = zeros(totalNumEns,1);
for ii = 1:totalNumEns
    tgSelector = cellTable.ensNum == ii & cellTable.tgCell == 1;
    
    numTg(ii) = sum(tgSelector);
    tgRespAve(ii) = nanmean(cellTable.dff(tgSelector));
    tgRespErr(ii) = nanstd(cellTable.dff(tgSelector))/sqrt(sum(tgSelector));
    
    numTgAct(ii) = sum(cellTable.dff(tgSelector)>0.5);
end

tgSelectorAll = cellTable.tgCell == 1;
tgRespAll = cellTable.dff(tgSelectorAll);

% Plot the results
figure(); clf;
subplot(3,2,1)
boxplot(tgRespAll)
ylabel('\DeltaF/F')
set(gca,'fontsize',16)

subplot(3,2,2)
histogram(tgRespAll)
xlabel('\DeltaF/F')
set(gca,'fontsize',16)

subplot(3,1,2); hold on;
hTemp = histogram(numTg);
hTemp2 = histogram(numTgAct);
plot(mean(numTg),max(hTemp.Values),'*','markersize',16,'color',[0 0.447 0.741])
plot(mean(numTgAct),max(hTemp2.Values),'*','markersize',16,'color',[0.85 0.325 0.098])
xlabel('# Matched Targets')
set(gca,'fontsize',16)
legend('All Tgs','dF/F>0.5 Tgs')

subplot(3,1,3)
errorbar(tgRespAve,tgRespErr)
set(gca,'fontsize',16)
ylabel('Ens Mean \DeltaF/F')
xlabel('Ens Number')

%% Plot ensemble response against ensemble response
ensResp = zeros(totalNumEns,1);
for ii = 1:totalNumEns
   cellSelector =  cellTable.ensNum == ii & cellTable.offTarget == 0 ...
       & cellTable.cellDist>50 & cellTable.cellDist<100;
   
   ensResp(ii) = nanmean(cellTable.dff(cellSelector));    
end

figure(); clf; hold on;
plot(tgRespAve,ensResp,'.','markersize',16)
plot([0 5], 0+[0 0],'k--','linewidth',1.5) 
set(gca,'fontsize',16)

fitLM2 = fit(tgRespAve,ensResp,'poly1');
plot(tgRespAve',fitLM2.p2+fitLM2.p1*tgRespAve','linewidth',2)
xlabel('Mean Ens \Delta F/F')
ylabel('Mean evoked \Delta F/F')

end

