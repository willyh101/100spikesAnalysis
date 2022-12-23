opts.distType = 'min';
opts.distBins =15:15:150; 
opts.distAxisRange = [0 150]; 

opts.plotTraces = 0;
opts.useVisCells = 0;
opts.useTunedCells = 0; %don't use tuned without vis

opts.variableCellFun = '';

figure(121); clf
lim =[-0.1 0.125]; [-0.4 0.25];

highMeanThresh = 0.5; %0.5; % 0.4687;
lowMeanThresh = 0.5;

lowEnsThresh = 0.3;
highEnsThresh = 0.7;

% Use these values for outVars.ensMaxD
% closeVal =  400;
% farVal = 500;

% Use these values for outVars.ensMeaD
closeVal = 200;
farVal = 200;

outInfo=[];
axs = [];
ax = subplot(2,2,4);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<lowEnsThresh & outVars.meanEnsOSI<lowMeanThresh;
% opts.criteriaToSplit =  outVars.ensMaxD;
opts.criteriaToSplit =  outVars.ensMeaD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('grey');
axs = [axs ax];

ax = subplot(2,2,2);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('grey');
axs = [axs ax];


% ax = subplot(3,2,3);
% opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI<lowEnsThresh & outVars.meanEnsOSI>highMeanThresh;
% % opts.criteriaToSplit = outVars.ensMaxD;
% opts.criteriaBins = [0 closeVal];% inf];
% [e outInfo{end+1}]= plotDistByCriteriaAx(All,outVars,opts,ax);
% e{1}.Color = rgb('sienna');
% axs = [axs ax];
% 
% ax = subplot(3,2,4);
% opts.criteriaBins = [farVal inf];% inf];
% [e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
% e{1}.Color = rgb('sienna');
% axs = [axs ax];


ax = subplot(2,2,3);
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10  &  outVars.ensOSI>highEnsThresh & outVars.meanEnsOSI>highMeanThresh;
% opts.criteriaToSplit = outVars.ensMaxD;
opts.criteriaBins = [0 closeVal];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('Magenta');
axs = [axs ax];

ax = subplot(2,2,1);
opts.criteriaBins = [farVal inf];% inf];
[e outInfo{end+1}] = plotDistByCriteriaAx(All,outVars,opts,ax);
e{1}.Color = rgb('Magenta');
axs = [axs ax];

linkaxes(axs)
ylim(lim);

disp('pVal first point diff from zero')
for i =1:4
    disp(num2str(signrank(outInfo{i}{1}.dat(:,1),0)))
end

[p h] = ranksum(outInfo{3}{1}.dat(:,1),outInfo{4}{1}.dat(:,1));
disp(['Tuned Near vs Far p= ' num2str(p)]);

[p h] = ranksum(outInfo{1}{1}.dat(:,1),outInfo{3}{1}.dat(:,1));
disp(['Near Untuned vs Tuned p= ' num2str(p)]);

[p h] = ranksum(outInfo{2}{1}.dat(:,1),outInfo{4}{1}.dat(:,1));
disp(['Far Untuned vs Tuned p= ' num2str(p)]);

[p h] = ranksum(outInfo{1}{1}.dat(:,1),outInfo{2}{1}.dat(:,1));
disp(['Untuned Near vs Far p= ' num2str(p)]);