opts.ensemblesToPlot = outVars.ensemblesToUse; 
opts.useVisCells =0;
opts.useTunedCells =0;
opts.minNumberOfCellsPerCondition = -1;

lowEnsThresh = 0.3;
highEnsThresh = 0.7;

highMeanThresh = 0.5;
lowMeanThresh = 0.5;

% Use these values for outVars.ensMeaD
closeVal = 200;
farVal = 200;

%% Close
opts.ensemblesToPlot = outVars.ensemblesToUse...
    & outVars.ensOSI>highEnsThresh...
    & outVars.meanEnsOSI>highMeanThresh...
    & outVars.ensMeaD < closeVal... 
    ; 

% iso
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & outVars.osi{ind}>0.25 & outVars.ensOriDiff{ens}==0';
[test2] = subsetPopResponse(All,outVars,opts);

closeIso = test2(opts.ensemblesToPlot);

% 45
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & outVars.osi{ind}>0.25 & outVars.ensOriDiff{ens}==45';
[test2] = subsetPopResponse(All,outVars,opts);

closeBtwn = test2(opts.ensemblesToPlot);

% Ortho
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & outVars.osi{ind}>0.25 & outVars.ensOriDiff{ens}==90';
[test2] = subsetPopResponse(All,outVars,opts);

closeOrtho = test2(opts.ensemblesToPlot);

% Non-Vis
opts.variableCellFun =  'outVars.pVisR{ind} > 0.05';
[test2] = subsetPopResponse(All,outVars,opts);

closeNonVis = test2(opts.ensemblesToPlot);

% Non-Tuned
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & (outVars.osi{ind}<=0.25 | isnan(outVars.ensOriDiff{ens}) )';
[test2] = subsetPopResponse(All,outVars,opts);

closeNonTuned = test2(opts.ensemblesToPlot);

closeData = [closeIso; closeBtwn; closeOrtho; closeNonVis; closeNonTuned];

figure(1);clf
title('Close')
plotSpread(closeData','showMM',4)
figure(2);clf
anova1(closeData')

%% Close
opts.ensemblesToPlot = outVars.ensemblesToUse...
    & outVars.ensOSI>highEnsThresh...
    & outVars.meanEnsOSI>highMeanThresh...
    & outVars.ensMeaD > farVal... 
    ; 

% iso
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & outVars.osi{ind}>0.25 & outVars.ensOriDiff{ens}==0';
[test2] = subsetPopResponse(All,outVars,opts);

farIso = test2(opts.ensemblesToPlot);

% 45
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & outVars.osi{ind}>0.25 & outVars.ensOriDiff{ens}==45';
[test2] = subsetPopResponse(All,outVars,opts);

farBtwn = test2(opts.ensemblesToPlot);

% Ortho
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & outVars.osi{ind}>0.25 & outVars.ensOriDiff{ens}==90';
[test2] = subsetPopResponse(All,outVars,opts);

farOrtho = test2(opts.ensemblesToPlot);

% Non-Vis
opts.variableCellFun =  'outVars.pVisR{ind} > 0.05';
[test2] = subsetPopResponse(All,outVars,opts);

farNonVis = test2(opts.ensemblesToPlot);

% Non-Tuned
opts.variableCellFun =  'outVars.pVisR{ind} < 0.05 & (outVars.osi{ind}<=0.25 | isnan(outVars.ensOriDiff{ens}) )';
[test2] = subsetPopResponse(All,outVars,opts);

farNonTuned = test2(opts.ensemblesToPlot);

farData = [farIso; farBtwn; farOrtho; farNonVis; farNonTuned];

figure(3);clf
title('far')
plotSpread(farData','showMM',4)
figure(4);clf
anova1(farData')