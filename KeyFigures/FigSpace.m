function [ensResp] = FigSpace(cellTable)

ensResp = zeros(160,1);
ensSpread = zeros(160,1);
for ii = 1:160
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

