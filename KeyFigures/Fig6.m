%%
% Reproduces close vs. far co-tuned plots
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig6(cellTable,cellCond)

% ensDistMetric = cellTable.cellEnsMaxD;
% spatialThresh = [-inf 400; 500 inf];
ensDistMetric = cellTable.cellEnsMeaD;
spatialThresh = [-inf 200; 200 inf];

totalNumEns = cellTable.ensNum(end);
distBins = [15:15:150];
plotDist = distBins(1:end-1) + 15/2;

% Ensemble thresholds
ensThreshs = [-inf 0.3; -inf 0.3;  0.7 inf];
meanEnsThreshs = [-inf 0.5; 0.5 inf; 0.5 inf];

% Conditions for this analysis
ensSelectorSpreadTight = ensDistMetric>spatialThresh(1,1) & ensDistMetric<spatialThresh(1,2);
ensSelectorSpreadLoose = ensDistMetric>spatialThresh(2,1) & ensDistMetric<spatialThresh(2,2);

respAveTight = zeros(length(distBins)-1,3); respAveLoose = zeros(length(distBins)-1,3);
respStdErrTight = zeros(length(distBins)-1,3); respStdErrLoose = zeros(length(distBins)-1,3);
numEnsUsedTight = zeros(length(distBins)-1,3); numEnsUsedLoose = zeros(length(distBins)-1,3);

% Preallocation
tightFirst = zeros(3,totalNumEns);
looseFirst = zeros(3,totalNumEns);

% Loop over conditions
for jj = 1:3
    ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(jj,1) & cellTable.cellEnsOSI<ensThreshs(jj,2)...
        & cellTable.cellMeanEnsOSI>meanEnsThreshs(jj,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(jj,2);
    
    cellDataAveTight=zeros(length(distBins)-1,totalNumEns);
    cellDataAveLoose=zeros(length(distBins)-1,totalNumEns);
    % Loop over ensembles
    for ii = 1:totalNumEns
        % Loop over distances
        for ll = 1:length(distBins)-1
            cellSelectorDist = cellTable.ensNum == ii & ...
                cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);
           
            cellSelectorTight = cellSelectorDist & ensSelectorTuning & ensSelectorSpreadTight & cellCond;
            cellSelectorLoose = cellSelectorDist & ensSelectorTuning & ensSelectorSpreadLoose & cellCond;
            
            cellDataAveTight(ll,ii) = nanmean(cellTable.dff(cellSelectorTight));
            cellDataAveLoose(ll,ii) = nanmean(cellTable.dff(cellSelectorLoose));
                        
            % Keep track of the number of ensembles used at each distance
            numEnsUsedTight(ll,jj) = numEnsUsedTight(ll,jj) + sign(sum(cellSelectorTight));
            numEnsUsedLoose(ll,jj) = numEnsUsedLoose(ll,jj) + sign(sum(cellSelectorLoose));
        end
    end
    
    % Averave across ensembles
    respAveTight(:,jj) = nanmean(cellDataAveTight,2);
    respStdErrTight(:,jj) = nanstd(cellDataAveTight,[],2)./sqrt(numEnsUsedTight(:,jj));
    
    respAveLoose(:,jj) = nanmean(cellDataAveLoose,2);
    respStdErrLoose(:,jj) = nanstd(cellDataAveLoose,[],2)./sqrt(numEnsUsedLoose(:,jj));
    
    %%
    tightFirst(jj,:) = cellDataAveTight(1,:);
    looseFirst(jj,:) = cellDataAveLoose(1,:);
end

%%
% strOps = {'','*'};
% if jj == 1
%     fprintf('1st point different from zero (untuned)\n')
% elseif jj == 2
%     fprintf('1st point different from zero (mixed tuned)\n')
% else
%     fprintf('1st point different from zero (co-tuned)\n')
% end
% p1 = signrank(cellDataAveTight(1,:),0); s1 = strOps{1+(p1<0.05)};
% p2 = signrank(cellDataAveLoose(1,:),0); s2 = strOps{1+(p2<0.05)};
% fprintf('Two-sided test, Tight: %.3f%s Loose: %.3f%s\n',p1,s1,p2,s2)
% 
% p2=signrank(cellDataAveLoose(1,:),0,'Tail','right'); s2 = strOps{1+(p2<0.05)};
% if jj < 3
%     p1=signrank(cellDataAveTight(1,:),0,'Tail','right'); s1 = strOps{1+(p1<0.05)};
% else
%     p1=signrank(cellDataAveTight(1,:),0,'Tail','left'); s1 = strOps{1+(p1<0.05)};
% end
% fprintf('One-sided test, Tight: %.3f%s Loose: %.3f%s\n',p1,s1,p2,s2)
% 
% p=ranksum(cellDataAveTight(1,:),cellDataAveLoose(1,:));
% fprintf('Tight vs Spread p=%.3f\n',p);
% 
% p=ranksum(cellDataAveTight(1,:),cellDataAveLoose(1,:),'tail','left');
% fprintf('Tight vs Spread p=%.3f\n',p);

%     [~,p]=ttest2(cellDataAveTight(1,:),cellDataAveLoose(1,:));

%%

% Plot the results
colorScheme =[];
colorScheme(1,1,:) = [0 0 0]; colorScheme(1,2,:) = [0 0 0]+0.5;
colorScheme(2,1,:) = [92, 64, 51]/255; colorScheme(2,2,:) = [165, 42, 42]/255;
colorScheme(3,1,:) = [0.494 0.184 0.556]; colorScheme(3,2,:) = [1 0 1];


%%


figure(3148866); clf; hold on;
bar([nanmedian(tightFirst(1,:)) nanmedian(tightFirst(3,:)) ...
    nanmedian(looseFirst(1,:)) nanmedian(looseFirst(3,:))])
% bar(nanmedian(tightFirst(1,:)),'facecolor',colorScheme(1,1,:))
plot(1,tightFirst(1,:),'.','markersize',16,'color',colorScheme(1,1,:))
plot(2,tightFirst(3,:),'.','markersize',16,'color',colorScheme(3,1,:))

plot(3,looseFirst(1,:),'.','markersize',16,'color',colorScheme(1,2,:))
plot(4,looseFirst(3,:),'.','markersize',16,'color',colorScheme(3,2,:))
set(gca,'fontsize',16)

%%
strOps = {'','*'};
fprintf('----------Rank sum test---------------\n')
p=ranksum(tightFirst(1,:),tightFirst(3,:)); s = strOps{1+(p<0.05)};
fprintf('Tight Untuned vs Tight co-tuned p=%.3f%s\n',p,s);
p=ranksum(looseFirst(1,:),looseFirst(3,:)); s = strOps{1+(p<0.05)};
fprintf('Loose Untuned vs Loose co-tuned p=%.3f%s\n',p,s);
p=ranksum(tightFirst(3,:),looseFirst(3,:)); s = strOps{1+(p<0.05)};
fprintf('Tight co-tuned vs Loose co-tuned p=%.3f%s\n',p,s);
p=ranksum(tightFirst(1,:),looseFirst(1,:)); s = strOps{1+(p<0.05)};
fprintf('Tight untuned vs Loose untuned p=%.3f%s\n',p,s);

%% t-test
% fprintf('----------T-test---------------\n')
% [~,p]=ttest2(tightFirst(1,:),tightFirst(3,:));
% fprintf('Tight Untuned vs Tight co-tuned p=%.3f\n',p);
% [~,p]=ttest2(looseFirst(1,:),looseFirst(3,:));
% fprintf('Loose Untuned vs Loose co-tuned p=%.3f\n',p);
% [~,p]=ttest2(tightFirst(3,:),looseFirst(3,:));
% fprintf('Tight co-tuned vs Loose co-tuned p=%.3f\n',p);
% [~,p]=ttest2(tightFirst(1,:),looseFirst(1,:));
% fprintf('Tight untuned vs Loose untuned p=%.3f\n',p);

%%
fprintf('-------------------------\n')
strOps = {'','*'};
p1 = signrank(tightFirst(1,:),0); s1 = strOps{1+(p1<0.05)};
p2 = signrank(tightFirst(3,:),0); s2 = strOps{1+(p2<0.05)};
p3 = signrank(looseFirst(1,:),0); s3 = strOps{1+(p3<0.05)};
p4 = signrank(looseFirst(3,:),0); s4 = strOps{1+(p4<0.05)};
fprintf('1st point different from zero (two-sided test) \n')
fprintf('Tight untuned: %.3f%s, Tight co-tuned: %.3f%s\n',p1,s1,p2,s2)
fprintf('Loose untuned: %.3f%s, Loose co-tuned: %.3f%s\n',p3,s3,p4,s4)

% plot(2,tightFirst(2,:),'k.','markersize',16)
% plot(3,tightFirst(3,:),'k.','markersize',16)

%%
figure(); clf; 
for jj = 1:2:3
    if jj == 1
        subplot(2,2,1); hold on;
    else
        subplot(2,2,3); hold on;
    end
    plot(plotDist,respAveTight(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,1,:))
    errorbar(plotDist,respAveTight(:,jj),respStdErrTight(:,jj),'linewidth',1.5,'color',colorScheme(jj,1,:),...
        'capsize',0)
    plot([0 250],0*[0 250],'k--')
    set(gca,'fontsize',16)
    xlim([0 150])
    ylim([-0.06 0.11])
    title(sprintf('Num Ens: %d',min(numEnsUsedTight(:,jj))))
    
    if jj == 1
        subplot(2,2,2); hold on;
    else
        subplot(2,2,4); hold on;
    end
    plot(plotDist,respAveLoose(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,2,:))
    errorbar(plotDist,respAveLoose(:,jj),respStdErrLoose(:,jj),'linewidth',1.5,'color',colorScheme(jj,2,:),...
        'capsize',0)
    plot([0 250],0*[0 250],'k--')
    set(gca,'fontsize',16)
    xlim([0 150])
    ylim([-0.06 0.11])
    title(sprintf('Num Ens: %d',min(numEnsUsedLoose(:,jj))))
end

%%
% figure(); clf; 
% for jj = 1:3
%     subplot(3,2,1+(jj-1)*2); hold on;
%     plot(plotDist,respAveTight(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,1,:))
%     errorbar(plotDist,respAveTight(:,jj),respStdErrTight(:,jj),'linewidth',1.5,'color',colorScheme(jj,1,:),...
%         'capsize',0)
%     plot([0 250],0*[0 250],'k--')
%     set(gca,'fontsize',16)
%     xlim([0 150])
%     ylim([-0.06 0.11])
%     title(sprintf('Num Ens: %d',min(numEnsUsedTight(:,jj))))
%     
%     subplot(3,2,2+(jj-1)*2); hold on;
%     plot(plotDist,respAveLoose(:,jj),'.-','linewidth',1.5,'markersize',15,'color',colorScheme(jj,2,:))
%     errorbar(plotDist,respAveLoose(:,jj),respStdErrLoose(:,jj),'linewidth',1.5,'color',colorScheme(jj,2,:),...
%         'capsize',0)
%     plot([0 250],0*[0 250],'k--')
%     set(gca,'fontsize',16)
%     xlim([0 150])
%     ylim([-0.06 0.11])
%     title(sprintf('Num Ens: %d',min(numEnsUsedLoose(:,jj))))
% end



%%


% disp('pVal first point diff from zero')
% for i =1:6
%     disp(num2str(signrank(outInfo{i}{1}.dat(:,1),0)))
% end
% 
% [p h] = ranksum(outInfo{5}{1}.dat(:,1),outInfo{6}{1}.dat(:,1));
% disp(['Tuned Near vs Far p= ' num2str(p)]);
% 
% [p h] = ranksum(outInfo{1}{1}.dat(:,1),outInfo{5}{1}.dat(:,1));
% disp(['Near Untuned vs Tuned p= ' num2str(p)]);
% 
% [p h] = ranksum(outInfo{2}{1}.dat(:,1),outInfo{6}{1}.dat(:,1));
% disp(['Far Untuned vs Tuned p= ' num2str(p)]);


end

