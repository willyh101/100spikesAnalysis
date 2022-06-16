%%
% Reproduces iso vs. ortho plot for close vs. far co-tuned condition
%
% cellCond is a vector of 1's and 0's that denotes which cells should be
% included (e.g., only non-offTarget cells)
%
% Should only included tuned and visually responsive cells
%
% Run cellByCellAnalysis_GH to use this function
%%
function Fig6IsoOrtho(cellTable,cellCondTuned,cellCondNonVis)


% ensDistMetric = cellTable.cellEnsMaxD;
% spatialThresh = [-inf 400; 500 inf];
ensDistMetric = cellTable.cellEnsMeaD;
spatialThresh = [-inf 200; 200 inf];

totalNumEns = cellTable.ensNum(end);
distBins = [15:15:150];
plotDist = distBins(1:end-1) + diff(distBins(1:2))/2;

% Ensemble thresholds
ensThreshs = [0.7 inf];
meanEnsThreshs = [0.5 inf];

% Conditions for this analysis
ensSelectorTuning = cellTable.cellEnsOSI>ensThreshs(1,1) & cellTable.cellEnsOSI<ensThreshs(1,2)...
    & cellTable.cellMeanEnsOSI>meanEnsThreshs(1,1) & cellTable.cellMeanEnsOSI<meanEnsThreshs(1,2);

cellSelectorOri = [cellTable.cellOrisDiff==0 cellTable.cellOrisDiff==45 ...
    cellTable.cellOrisDiff==90 cellTable.cellOrisDiff>0];

% Iso, 45, Ortho, Not Iso, Non-Vis
num_conds = 5;
respAve = zeros(length(distBins)-1,2,num_conds); 
respStdErr = zeros(length(distBins)-1,2,num_conds); 
numEnsUsed = zeros(length(distBins)-1,2,num_conds); 

% Loop over conditions
for jj = 1:2
    
    ensSelectorSpread = ensDistMetric>spatialThresh(jj,1) & ensDistMetric<spatialThresh(jj,2);
    
    
    cellDataAve = zeros(length(distBins)-1,totalNumEns,num_conds);
    % Loop over ensembles
    for ii = 1:totalNumEns
        % Loop over distances
        for ll = 1:length(distBins)-1
            cellSelectorDist = cellTable.ensNum == ii & ...
                cellTable.cellDist>distBins(ll) & cellTable.cellDist<distBins(ll+1);

            for gg = 1:5
                
                if gg < 5
                    cellSelector = cellSelectorDist & ensSelectorTuning & ensSelectorSpread ...
                        & cellSelectorOri(:,gg) & cellCondTuned;
                else
                    cellSelector = cellSelectorDist & ensSelectorTuning & ensSelectorSpread ...
                        & cellCondNonVis;
                end

                cellDataAve(ll,ii,gg) = nanmean(cellTable.dff(cellSelector));
                % Keep track of the number of ensembles used at each distance
                numEnsUsed(ll,jj,gg) = numEnsUsed(ll,jj,gg)+sign(sum(cellSelector));
            end
        end
    end
    
    % Averave across ensembles
    for gg = 1:5
        respAve(:,jj,gg) = nanmean(cellDataAve(:,:,gg),2);
        respStdErr(:,jj,gg) = nanstd(cellDataAve(:,:,gg),[],2)./sqrt(numEnsUsed(:,jj,gg));
        
        % Save first bin
        resp{jj,gg} = cellDataAve(1,~isnan(cellDataAve(1,:,gg)),gg);
    end    
end

%%
      
% [~,p]=ttest2(respIso{1},respIso{2})
% [~,p]=ttest2(resp45{1},resp45{2},'tail','left')
% [~,p]=ttest2(respOrtho{1},respOrtho{2},'tail','left')

%%

% [~,p] = ttest2(respIso{1},respOrtho{1},'tail','right')
% [~,p] = ttest2(respIso{2},respOrtho{2},'tail','right')
% 
% %%
% 
% [~,p]=ttest2(respIso{1},respIso{2})
% [~,p]=ttest2(respNotIso{1},respNotIso{2},'tail','left')
% %%
% [~,p]=ttest2(respIso{1},respNotIso{1},'tail','right')
% [~,p]=ttest2(respIso{2},respNotIso{2},'tail','right')
% 
% 
% [~,p]=ttest(respNotIso{2},0,'tail','left')
% [~,p]=ttest(respNotIso{1},0,'tail','left')
%% Plot the results
colorScheme =[];
colorScheme(1,:) = [0 0.447 0.741];
colorScheme(2,:) = [0.4940 0.184 0.556];
colorScheme(3,:) = [255, 0, 255]/255;
colorScheme(4,:) = [0.5 0.5 0.5];
titles={'Iso','45','Ortho','Non-Vis'};
indicesToPlot = [1 2 3 5];
numToPlot = length(indicesToPlot);

figure(668877); clf; 
for jj = 1:2
    for gg = 1:numToPlot
        subplot(2,numToPlot,gg+numToPlot*(jj-1)); hold on;
        errorbar(plotDist,respAve(:,jj,indicesToPlot(gg)),respStdErr(:,jj,indicesToPlot(gg)),...
            'linewidth',1.5,'color',colorScheme(gg,:),'capsize',0);
        xlim([0 150])
        plot([0 250],0*[0 250],'k--')
        ylim([-0.15 0.15])
        title(titles{gg})
        set(gca,'fontsize',16)
        if jj == 1 && gg == 1
            ylabel('Close')
        elseif jj == 2  && gg == 1
            ylabel('Far')
        end
    end
end

%%

figure(6688); clf; 

subplot(1,2,1); hold on;
% plot([respAveIso(1,1) respAve45(1,1) respAveOrtho(1,1)],'.','markersize',16,'color',[0    0.4470    0.7410])

errorbar(1:3,squeeze(respAve(1,1,1:3)),squeeze(respStdErr(1,1,1:3)),...
    'linewidth',1.5,'color',[0    0.4470    0.7410],'capsize',0)
errorbar(1:3,squeeze(respAve(1,2,1:3)),squeeze(respStdErr(1,2,1:3)),...
    'linewidth',1.5,'color',[0.8500    0.3250    0.0980],'capsize',0)
plot([1 4],0+[0 0],'k--')

errorbar(4,squeeze(respAve(1,1,5)),squeeze(respStdErr(1,1,5)),...
    'linewidth',1.5,'color',[0    0.4470    0.7410],'capsize',0)
errorbar(4,squeeze(respAve(1,2,5)),squeeze(respStdErr(1,2,5)),...
    'linewidth',1.5,'color',[0.8500    0.3250    0.0980],'capsize',0)

set(gca,'fontsize',16)
xticks([1 2 3 4])
xticklabels({'Iso', '', 'Ortho','Non-vis'})
xlabel('\Delta\theta')
ylabel('Mean Evoked \DeltaF/F')
xlim([0.75 4.25])
legend('Far','Close')


subplot(1,2,2); hold on;
diffAve = squeeze(respAve(1,1,[1 2 3 5]) -respAve(1,2,[1 2 3 5]));
diffError = sqrt(squeeze(respStdErr(1,1,[1 2 3 5])).^2+...
    squeeze(respStdErr(1,2,[1 2 3 5])).^2);
errorbar([1 2 3 4],diffAve,diffError,'linewidth',1.5,'capsize',0)
xlim([0.75 4.25])

plot([1 3],0+[0 0],'k--')
xticks([1 2 3 4])
xticklabels({'Iso', '', 'Ortho','Non-vis'})
set(gca,'fontsize',16)
ylabel('Far-Close \DeltaF/F')





%%
end

