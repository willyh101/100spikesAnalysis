%%
% Reproduces the effect of ensemble spread figure
%
% Run cellByCellAnalysis_GH to use this function
%%
function [ensResp] = FigSpace(cellTable)

totalNumEns = cellTable.ensNum(end);

ensResp = zeros(totalNumEns,1);
ensSpread = zeros(totalNumEns,1);
for ii = 1:totalNumEns
   cellSelector =  cellTable.ensNum == ii & cellTable.offTarget == 0 ...
       & cellTable.cellDist>50 & cellTable.cellDist<100;
   
   ensResp(ii) = nanmean(cellTable.dff(cellSelector));
   ensSpread(ii) = mean(cellTable.cellEnsMaxD(cellSelector));
    
end

figure(); hold on;
plot(ensSpread,ensResp,'.','markersize',16)
plot([150 650], 0+[0 0],'k--','linewidth',1.5) 
set(gca,'fontsize',16)
fitLM = fit(double(ensSpread),ensResp,'poly1');
plot(ensSpread,fitLM.p2+fitLM.p1*ensSpread,'linewidth',2)
xlabel('Ens Spread')
ylabel('Mean evoked \Delta F/F')

end

