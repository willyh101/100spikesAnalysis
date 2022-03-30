
muPerPx = 800/512;
opts.thisPlaneTolerance =15/muPerPx;% 15/muPerPx;
opts.onePlaneTolerance = 30/muPerPx; %30/muPerPx;

recalcOffTargetRiskNoArt;

%%s
opts.nearbyActiveZThreshold = 3;
opts.nearbyActiveThreshold =1; -inf; 0.5;



nearbySizes=[];
for i = 1:numel(outVars.offTargetRiskEns)
    otr = outVars.offTargetRiskEns{i};
    ind = outVars.ensIndNumber(i);
    
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    s = outVars.ensHNumber(i); %i think this is the order of resp
    h = All(ind).out.exp.stimParams.roi{s};
     htg = All(ind).out.exp.holoTargets{h};
    htg(isnan(htg))=[];
    
    meanVals = squeeze(respMat(s,1,:) - baseMat(s,1,:));
    
    nearbyList = find(zscore(meanVals)>opts.nearbyActiveZThreshold & otr');
            nearbyList = find(meanVals>opts.nearbyActiveThreshold & otr');
% 
    nearbyEns{i} = nearbyList;
    nearbySizes(i)  = numel(nearbyList);
    
    meanOSI(i) = nanmean(All(ind).out.anal.osi(nearbyList));
    
    % get mean Ori by simple mean
    meanOri(i) = nanmean(idx2ori(All(ind).out.anal.prefOri(nearbyList), [nan 0:45:315]));
    
    curve = nanmean(All(ind).out.anal.oriCurve(:, nearbyList), 2); %set in getTuningCurve
    [~, ensPOidx] = max(curve);
    curve(:, ensPOidx==1) = nan;
    curve(1, :) = [];
    ensPO(i) = idx2ori(ensPOidx, [nan 0:45:315]);
    ensOSI(i) = osi(curve);
end

outVars.ensNearbyOSI = ensOSI;
outVars.ensNearbyPO = ensPO;
outVars.meanNearbyOSI = meanOSI;
outVars.meanNearbyOri = meanOri;

figure(3);clf
subplot(1,3,1)
histogram(nearbySizes(outVars.ensemblesToUse))
title('Nearby ''ensemble'' sizes')
subplot(1,3,2)

plot(outVars.ensNearbyOSI(outVars.ensemblesToUse), outVars.ensOSI(outVars.ensemblesToUse),'o')
subplot(1,3,3)
histogram(ensOSI(outVars.ensemblesToUse))
title('Nearby ''ensemble'' ensOSI')


%%
muPerPx = 800/512;
opts.thisPlaneTolerance =15/muPerPx;% 15/muPerPx;
opts.onePlaneTolerance = 30/muPerPx; %30/muPerPx;

recalcOffTargetRisk;
opts.distType = 'min';
opts.distBins =15:10:150; %15:15:150; %[15:15:150]; [opts.thisPlaneTolerance*muPerPx:10:150]; %[0:25:150]; 
opts.plotTraces = 0;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
figure(19); clf
lim =[-0.1 0.2]; [-0.4 0.25];

meanThresh = 0.5; %0.5; % 0.4687;
closeVal =  400
farVal =500

outInfo=[];
axs = [];
ax = subplot(3,2,1);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensNearbyOSI<0.3 & outVars.meanNearbyOSI<meanThresh;
opts.criteriaToSplit =  outVars.ensMaxD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('grey');
axs = [axs ax];

ax = subplot(3,2,2);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('grey');
axs = [axs ax];


ax = subplot(3,2,3);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensNearbyOSI<0.3 & outVars.meanNearbyOSI>meanThresh;
% opts.criteriaToSplit = outVars.ensMaxD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('sienna');
axs = [axs ax];

ax = subplot(3,2,4);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('sienna');
axs = [axs ax];


ax = subplot(3,2,5);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensNearbyOSI>0.7 & outVars.meanNearbyOSI>meanThresh;
% opts.criteriaToSplit = outVars.ensMaxD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('Magenta');
axs = [axs ax];

ax = subplot(3,2,6);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('Magenta');
axs = [axs ax];

linkaxes(axs)
ylim(lim);

disp('pVal first point diff from zero')
for i =1:6
    disp(num2str(signrank(outInfo{i}{1}.dat(:,1),0)))
end

[p h] = ranksum(outInfo{5}{1}.dat(:,1),outInfo{6}{1}.dat(:,1));
disp(['Tuned Near vs Far p= ' num2str(p)]);

[p h] = ranksum(outInfo{1}{1}.dat(:,1),outInfo{5}{1}.dat(:,1));
disp(['Near Untuned vs Tuned p= ' num2str(p)]);

[p h] = ranksum(outInfo{2}{1}.dat(:,1),outInfo{6}{1}.dat(:,1));
disp(['Far Untuned vs Tuned p= ' num2str(p)]);