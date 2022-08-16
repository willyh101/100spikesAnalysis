%% Spread vs effect by distance bin 
%%Close Excitation
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;

opts.variableCellFun =  '(outVars.distToEnsemble{i}<30)';
% opts.variableCellFun =  '(outVars.distToEnsemble{i}>50 & outVars.distToEnsemble{i}<100)';
[midResponse] = subsetPopResponse(All,outVars,opts);

figure(19);clf
s1 = subplot(1,3,1);

% x = outVars.ensOSI(opts.ensemblesToPlot)';  
x = outVars.ensMeaD(opts.ensemblesToPlot)';  
y = midResponse(opts.ensemblesToPlot);

scatter(x,y,[],rgb('LimeGreen'),'filled')
refline(0)
xlabel('Spread ') ;
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');
hold on
plot(fs)
legend off

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp('Close: ')
disp(['slope is: ' num2str(fs.p1)])
disp(['Pval is: ' num2str(pVal(1))])

%%Mid Distance Inhibition
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;

% opts.variableCellFun =  '(outVars.distToEnsemble{i}<30)';
opts.variableCellFun =  '(outVars.distToEnsemble{i}>50 & outVars.distToEnsemble{i}<100)';
[midResponse] = subsetPopResponse(All,outVars,opts);

figure(19);
s2 = subplot(1,3,2);

% x = outVars.ensOSI(opts.ensemblesToPlot)';  
x = outVars.ensMeaD(opts.ensemblesToPlot)';  
y = midResponse(opts.ensemblesToPlot);

scatter(x,y,[],rgb('SteelBlue'),'filled')
refline(0)
xlabel('Spread ') ;
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');
hold on
plot(fs)
legend off

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp('Mid: ')
disp(['slope is: ' num2str(fs.p1)])
disp(['Pval is: ' num2str(pVal(1))])

%%Far Distance Inhibition
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;

% opts.variableCellFun =  '(outVars.distToEnsemble{i}<30)';
opts.variableCellFun =  '(outVars.distToEnsemble{i}>150)';
[midResponse] = subsetPopResponse(All,outVars,opts);

figure(19);
s3 = subplot(1,3,3);


% x = outVars.ensOSI(opts.ensemblesToPlot)';  
x = outVars.ensMeaD(opts.ensemblesToPlot)';  
y = midResponse(opts.ensemblesToPlot);

scatter(x,y,[],rgb('Amethyst'),'filled')
refline(0)
xlabel('Spread ') ;
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');
hold on
plot(fs)
legend off

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp('Close: ')
disp(['slope is: ' num2str(fs.p1)])
disp(['Pval is: ' num2str(pVal(1))])
%%meta
linkaxes([s1,s2,s3], 'x');
xlim([95 325])
s2.YLim=[-0.08 0.08];
s3.YLim=[-0.08 0.08];

%% Now based on ensOSI
%%Close Excitation
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;

opts.variableCellFun =  '(outVars.distToEnsemble{i}<30)';
% opts.variableCellFun =  '(outVars.distToEnsemble{i}>50 & outVars.distToEnsemble{i}<100)';
[midResponse] = subsetPopResponse(All,outVars,opts);

figure(19);clf
s1 = subplot(1,3,1);

x = outVars.ensOSI(opts.ensemblesToPlot)';  
% x = outVars.ensMeaD(opts.ensemblesToPlot)';  
y = midResponse(opts.ensemblesToPlot);

scatter(x,y,[],rgb('LimeGreen'),'filled')
refline(0)
xlabel('Spread ') ;
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');
hold on
plot(fs)
legend off

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp('Close: ')
disp(['slope is: ' num2str(fs.p1)])
disp(['Pval is: ' num2str(pVal(1))])

%%Mid Distance Inhibition
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;

% opts.variableCellFun =  '(outVars.distToEnsemble{i}<30)';
opts.variableCellFun =  '(outVars.distToEnsemble{i}>50 & outVars.distToEnsemble{i}<100)';
[midResponse] = subsetPopResponse(All,outVars,opts);

figure(19);
s2 = subplot(1,3,2);

x = outVars.ensOSI(opts.ensemblesToPlot)';  
% x = outVars.ensMeaD(opts.ensemblesToPlot)';  
y = midResponse(opts.ensemblesToPlot);

scatter(x,y,[],rgb('SteelBlue'),'filled')
refline(0)
xlabel('Spread ') ;
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');
hold on
plot(fs)
legend off

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp('Mid: ')
disp(['slope is: ' num2str(fs.p1)])
disp(['Pval is: ' num2str(pVal(1))])

%%Far Distance Inhibition
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;

% opts.variableCellFun =  '(outVars.distToEnsemble{i}<30)';
opts.variableCellFun =  '(outVars.distToEnsemble{i}>150)';
[midResponse] = subsetPopResponse(All,outVars,opts);

figure(19);
s3 = subplot(1,3,3);


x = outVars.ensOSI(opts.ensemblesToPlot)';  
% x = outVars.ensMeaD(opts.ensemblesToPlot)';  
y = midResponse(opts.ensemblesToPlot);

scatter(x,y,[],rgb('Amethyst'),'filled')
refline(0)
xlabel('Spread ') ;
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');
hold on
plot(fs)
legend off

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp('Close: ')
disp(['slope is: ' num2str(fs.p1)])
disp(['Pval is: ' num2str(pVal(1))])
%%meta
linkaxes([s1,s2,s3], 'x');
% xlim([95 325])
s2.YLim=[-0.08 0.08];
s3.YLim=[-0.08 0.08];

%% all Data
opts.ensemblesToPlot = outVars.ensemblesToUse & outVars.numCellsEachEnsBackup==10;
opts.useVisCells = 0;
opts.useTunedCells =0; %don't use tuned without vis
opts.minNumberOfCellsPerCondition = -1;

% opts.variableCellFun =  '(outVars.distToEnsemble{i}<30)';
opts.variableCellFun =  '(outVars.distToEnsemble{i}>15)';
[midResponse] = subsetPopResponse(All,outVars,opts);

figure(20);
% s3 = subplot(1,3,3);


x = outVars.ensOSI(opts.ensemblesToPlot)';  
% x = outVars.ensMeaD(opts.ensemblesToPlot)';  
y = midResponse(opts.ensemblesToPlot);

scatter(x,y,[],rgb('black'),'filled')
refline(0)
xlabel('Spread ') ;
nanEither = isnan(x) | isnan(y');

[fs, gs] = fit(x(~nanEither),y(~nanEither)','poly1');
hold on
plot(fs)
legend off

[p Rsq pVal] = simplifiedLinearRegression(x(~nanEither),y(~nanEither)');
disp('Close: ')
disp(['slope is: ' num2str(fs.p1)])
disp(['Pval is: ' num2str(pVal(1))])