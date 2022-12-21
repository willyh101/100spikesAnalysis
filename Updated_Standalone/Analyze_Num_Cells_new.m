%%
% Same as Analyze_Feature except with the new pre-processing step of
% removing offTarget cells from the previous trial from analysis
%
% New code block starts at line ~71
% Also change meanMatrixVisandCorr to meanMatrixVisandCorr_v2 (replaces
% some mean() functions with nanmean()
%%
clear; close all; clc;

%% loadLists

% oriLoadList;
allLoadList;

% loadPath = 'path/to/outfiles/directory';
% loadPath = 'T:\Outfiles';
loadPath = '/Users/greghandy/Research_Local_v2/'; % where ever your files are


addpath(genpath(loadPath))

%% Load data

numExps = numel(loadList);
disp(['There are ' num2str(numExps) ' Exps in this LoadList'])
if numExps ~= 0
    clear All
    if ~iscell(loadList)
        numExps=1;
        temp = loadList;
        clear loadList;
        loadList{1} = temp;
    end
    for ind = 1:numExps
        pTime =tic;
        fprintf(['Loading Experiment ' num2str(ind) '...']);
        All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
        fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
    end
else
    disp('Did you press this by accident?')
end

%% error fixer
%CAUTION! ERRORS WILL OCCUR IF YOU RUN MORE THAN ONCE!
[All] = allLoadListErrorFixer(All,loadList);

%% Rematch Targeted Cells
opts.matchDistanceMicrons = 12; 15.6; %6;
rematchTargetedCells;

%% clean Data, and create fields.
opts.FRDefault=6;
opts.recWinRange = [0.5 1.5]; [0.5 1.5];%[1.5 2.5];%[0.5 1.5];% %from vis Start in s [1.25 2.5];

%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull
opts.stimSuccessByZ = 1; %if 0 calculates based on dataTouse, if 1 calculates based on zdfDat;

%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded
opts.runValPercent = 0.75; %percent of frames that need to be below run threshold

%% New Code for the removal of offTarget cells 

% Note: cleanData is run here and after dataToUse code block
[All, ~] = cleanData(All,opts);

% Calculate the offTarget risk
muPerPx=800/512;
opts.thisPlaneTolerance = 15/muPerPx; %11.25;%7.5;%1FWHM%10; %in pixels
opts.onePlaneTolerance = 30/muPerPx; %22.5;%15;%2FWHM %20;
[All] = calcOffTargetRisk(All,opts);

% Remove cells from analysis if they were an offTarget in the 
% previous trial
[All] = prevOffTargetRemoval(All);

%% Set Data To use
for ind=1:numExps
    All(ind).out.exp.dataToUse = All(ind).out.exp.dfData;
end
disp('Data To Use is set')

%% Re-run clean data, and create fields.
[All, outVars] = cleanData(All,opts);

ensStimScore        = outVars.ensStimScore;
hzEachEns           = outVars.hzEachEns;
numCellsEachEns     = outVars.numCellsEachEns;
numSpikesEachStim   = outVars.numSpikesEachStim;
percentLowRunTrials = outVars.percentLowRunTrials;
numSpikesEachEns    = outVars.numSpikesEachEns;
numSpikesEachCell   = outVars.numSpikesEachCell;

outVars.numCellsEachEnsBackup = outVars.numCellsEachEns;

names=[];
for Ind = 1:numel(All)
    names{Ind}=lower(strrep(All(Ind).out.info.mouse, '_', '.'));
end
outVars.names = names;

%% restrict Cells to use
opts.minMeanThreshold = 0.25;
opts.maxMeanThreshold = inf;

opts.verbose =0;
[All, cellExcludeResults] = cellExcluder(All,opts); 
allResults = cat(1,cellExcludeResults{:});
disp(['In total ' num2str(sum(allResults)) ' Cells Excluded. ' num2str(mean(allResults)*100,2) '%']);
disp(['Overall ' num2str(sum(~allResults)) ' Cells Passed!'])

opts.minNumCellsInd=250;
tooFewCellsInds = cellfun(@(x) sum(~x)<opts.minNumCellsInd,cellExcludeResults);
disp([ num2str(sum(tooFewCellsInds)) ' inds have < ' num2str(opts.minNumCellsInd) ' cells, and should be exccluded']);



%% Make all dataPlots into matrixes of mean responses
%%Determine Vis Responsive and Process Correlation

opts.visAlpha = 0.05;

%oftarget risk params
opts.thisPlaneTolerance = 11.25;%7.5;%1FWHM%10; %in um;% pixels
opts.onePlaneTolerance = 22.5;%15;%2FWHM %20;
opts.distBins =  [0:25:1000]; [0:25:1000];
opts.skipVis =1;

[All, outVars] = meanMatrixVisandCorr_v2(All,opts,outVars); %one of the main analysis functions

visPercent = outVars.visPercent;
outVars.visPercentFromExp = visPercent;
ensIndNumber =outVars.ensIndNumber;



%% Misc Additional Variables:
%%RedSection: if there is a red section (will run even if not...)
[outVars] = detectShotRedCells(All,outVars);
ensHasRed = outVars.ensHasRed;

try
arrayfun(@(x) sum(~isnan(x.out.red.RedCells)),All)
arrayfun(@(x) mean(~isnan(x.out.red.RedCells)),All)
catch end
%%Identify the Experiment type for comparison or exclusion
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;
ensExpressionType = indExpressionType(outVars.ensIndNumber);
outVars.ensExpressionType = ensExpressionType;

%%Missed Target Exclusion Criteria
%detects if too many targets were not detected in S2p
opts.FractionMissable = 0.33; %what percent of targets are missable before you exclude the ens
[outVars] = missedTargetDetector(All,outVars,opts);

ensMissedTargetF = outVars.ensMissedTargetF; %Fraction of targets per ensemble Missed
ensMissedTarget = outVars.ensMissedTarget; %Ensemble is unuseable
numMatchedTargets = outVars.numMatchedTargets;
%%Determine Date
ensDate=[];
for i = 1:numel(ensIndNumber)
    ensDate(i) = str2num(All(ensIndNumber(i)).out.info.date);
end
outVars.ensDate=ensDate;

%%Identify duplicate holograms
[outVars] = identifyDuplicateHolos(All,outVars);

%% main Ensembles to Use section

[ensemblesToUse, outVars, opts] = ensemblesToUse_fn(All,numExps,...
    numSpikesEachStim,ensIndNumber,visPercent,percentLowRunTrials,...
    tooFewCellsInds,outVars,opts,ensExpressionType,numSpikesEachEns,...
    ensStimScore,ensHasRed,ensMissedTarget,numMatchedTargets,...
    numCellsEachEns,ensDate,names,'numCells');

%% Ensure Ensembles come from inds where each type of holos was presented.

ensSizes =[3 5 10 20 33]; [3 10 33];[5 10 20];
indsBySize =[];
for i = 1:numel(ensSizes);
    indsBySiz{i} = unique(ensIndNumber(ensemblesToUse & numCellsEachEns ==ensSizes(i)));
end
AllList = 1:numel(All);
temp = cellfun(@(x) ismember(AllList,x),indsBySiz,'uniformoutput',0);
indsWithAllSize = find(all(cat(1,temp{:})));

indsWithAllSize = find(sum(cat(1,temp{:}))>=2);

ensemblesToUse = ensemblesToUse & ismember(ensIndNumber,indsWithAllSize);
disp(['Total of ' num2str(sum(ensemblesToUse)) ' Ensembles Included.'])
disp([ num2str(numel(unique(ensIndNumber(ensemblesToUse)))) ' FOVs'])
disp([ num2str(numel(unique(names(unique(ensIndNumber(ensemblesToUse)))))) ' Mice']);

outVars.ensemblesToUse = ensemblesToUse;
outVars.IndsUsed = unique(ensIndNumber(ensemblesToUse));

%% Basic Response Plots
outVars.defaultColorMap = 'viridis';
plotResponseBySize(outVars,0)
ylim([-0.1 0.05])
yticks([-0.1:0.05:0.05])
