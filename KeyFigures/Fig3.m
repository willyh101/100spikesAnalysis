%%
% Reproduces the effect of ensemble spread figure
%
% Run expFigs2Through5_Dec22 to use this function
%%
function [ensResp] = Fig3(cellTable,cellCond)

ensDistMetric = cellTable.cellEnsMeaD;

totalNumEns = cellTable.ensNum(end);

ensResp = zeros(totalNumEns,1);
ensSpread = zeros(totalNumEns,1);
for ii = 1:totalNumEns
   cellSelector =  cellTable.ensNum == ii ...
       & cellTable.cellDist>50 & cellTable.cellDist<150 & cellCond;
   ensResp(ii) = nanmean(cellTable.dff(cellSelector));
   ensSpread(ii) = unique(ensDistMetric(cellSelector));
end

figure(); hold on;
plot(ensSpread,ensResp,'.','markersize',16)
plot([90 350], 0+[0 0],'k--','linewidth',1.5) 
set(gca,'fontsize',16)
fitLM = fit(double(ensSpread),ensResp,'poly1');
plot([90 350],fitLM.p2+fitLM.p1*[90 350],'linewidth',2)
xlabel('Ens Spread')
ylabel('Mean evoked \Delta F/F')
xlim([90 350])

%%
fprintf('Slope %e\n',fitLM.p1) 

%%
[~,~, pVal] = simplifiedLinearRegression(double(ensSpread),ensResp);
fprintf('p-val: %e \n',pVal(1))

end

