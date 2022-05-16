%% Load Experiments/Setup
clear
close all

addpath(genpath('100spikesAnalysis'), genpath('Ian Code'), genpath('analysis-code/matlab')) % will pathing
%%
% loadList = loadList(15);


% allLoadList;
% oriLoadList;
% SSTOriLoadList;
% PVOriLoadList;
% u19LoadList;
% manifoldLoadList;
ai203LoadList;

% loadPath = 'U:\ioldenburg\outputdata1'
% loadPath = 'C:\Users\ian\Dropbox\Adesnik\Data\outputdata1'
% loadPath = 'C:\Users\SabatiniLab\Dropbox\Adesnik\Data\outputdata1' %Ian Desktop
% loadPath = 'C:\Users\Will\Local Data\100spikes-results\outfiles-ori'
% loadPath = 'T:\Outfiles';
% loadPath = 'E:\100spikes-results\outfiles-master';
loadPath = 'E:\outfiles';
% loadPath = '/Volumes/Frankenshare/Outfiles';
% loadPath = '/Users/willh/Desktop/ai203_data';
addpath(genpath(loadPath))

%%

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

%% Fix Known Errors
%CAUTION! ERRORS WILL OCCUR IF YOU RUN MORE THAN ONCE!
[All] = allLoadListErrorFixer(All,loadList);

%% clean Data, and create fields.

opts.FRDefault=6;
opts.recWinRange = [0.5 1.5];% %from vis Start in s [1.25 2.5];


%Stim Success Thresholds
opts.stimsuccessZ = 0.25; %0.3, 0.25 over this number is a succesfull stim
opts.stimEnsSuccess = 0.5; %0.5, fraction of ensemble that needs to be succsfull

%run Threshold
opts.runThreshold = 6 ; %trials with runspeed below this will be excluded


% All = cleanDataVisOnly(All,opts);
All = cleanDataGMNOnly(All,opts); % assumes that GMN is vis2

names=[];
for Ind = 1:numel(All)
    names{Ind}=lower(strrep(All(Ind).out.info.mouse, '_', '.'));
end
outVars.names = names;

%% Make all dataPlots into matrixes of mean responses
%%Determine Vis Responsive and Process Correlation

opts.visAlpha = 0.05;

% opts.distBins =  [0:25:1000]; [0:25:1000];

% [All, outVars] = meanMatrixVisandCorr(All,opts,outVars); %one of the main analysis functions

% visPercent = outVars.visPercent;
% outVars.visPercentFromExp = visPercent;
% ensIndNumber =outVars.ensIndNumber;

%% vis
opts.visAlpha = 0.05;
% [All, outVars] = CalcPVisRFromVis(All,opts,outVars);
[All, outVars] = CalcPVisRFromGMN(All,opts,outVars); % again assumes that GMN is vis2
visPercent = outVars.visPercent;
outVars.visPercentFromVis = visPercent;

%% Identify the Experiment type for comparison or exclusion
[All,outVars] = ExpressionTypeIdentifier(All,outVars);
indExpressionType = outVars.indExpressionType;

%% CRFs

clear CRF CRF203
k = 0;
kk = 0;

for ind=1:numExps
    clear crf
    c = 0;
    if isfield(All(ind).out, 'vis2')
        c=c+1;
        
        vs = unique(All(ind).out.vis2.visID);
        pVisR = All(ind).out.anal.pVisR;
        winToUse = All(ind).out.vis2.win;
        bwinToUse = All(ind).out.vis2.bwin;

        All(ind).out.info
        pause

        for i=1:numel(vs)
            v = vs(i);

            cellsToUse = pVisR < 0.05;
            trialsToUse = All(ind).out.vis2.visID == v;

            if size(trialsToUse,2) > size(All(ind).out.vis2.rdata,2)
                trialsToUse = trialsToUse(1:size(All(ind).out.vis2.rdata,2));
            end

            resp = squeeze(mean(mean(All(ind).out.vis2.zdfData(cellsToUse, winToUse(1):winToUse(2), trialsToUse),2),3));
            base = squeeze(mean(mean(All(ind).out.vis2.zdfData(cellsToUse, bwinToUse(1):bwinToUse(2), trialsToUse),2),3));

            crf(i,:) = resp - base;
        end

        All(ind).out.anal.CRF = crf;

        if strcmp(All(ind).out.info.ExpressionType, 'Ai203')
            k = k+1;
            CRF203{k} = crf;
        else
            kk = kk+1;
            CRF{kk} = crf;
        end
    end 
end

crf = cell2mat(CRF);
crf2 = cell2mat(CRF203);
%% Plot CRFs by  cell

figure(1)
clf
hold on

% regular
% m = mean(crf,2);
% m2 = mean(crf2, 2);

% normalized
m = mean(crf,2);
m = m/max(m);
m2 = mean(crf2, 2);
m2 = m2/max(m2);

e = std(crf,[],2)./sqrt(size(crf,2));
e2 = std(crf2,[],2)./sqrt(size(crf2,2));

xs = [0, 1, 3, 10, 30, 100];

e = errorbar(xs, m, e);
e.LineWidth = 2;
e = errorbar(xs, m2, e2);
e.LineWidth = 2;

legend('tetO', 'Ai203')
% ylabel('Normalized Resp.')
ylabel('ZDF')
xlabel('% Contrast')
% set(gca, 'XSCale', 'log')
title('CRF by Cell')

n=size(crf,2);
n2=size(crf2,2);

disp(['tetO:  n = ' num2str(n) ' by cell'])
disp(['Ai203: n = ' num2str(n2) ' by cell'])


%% Plot CRFs  by FOV

figure(2)
clf
hold  on

dat = cell2mat(cellfun(@(x) mean(x,2), CRF, 'UniformOutput',0));
dat2 = cell2mat(cellfun(@(x) mean(x,2), CRF203, 'UniformOutput',0));

% regular
% m = mean(dat,2);
% m2 = mean(dat2, 2);

% normalized
m = mean(dat,2);
m = m/max(m);
m2 = mean(dat2, 2);
m2 = m2/max(m2);

e = std(dat,[],2)./sqrt(size(dat,2));
e2 = std(dat2,[],2)./sqrt(size(dat2,2));

xs = [0, 1, 3, 10, 30, 100];

e = errorbar(xs, m, e);
e.LineWidth = 2;
e = errorbar(xs, m2, e2);
e.LineWidth = 2;

legend('tetO', 'Ai203')
% ylabel('Normalized Resp.')
ylabel('ZDF')
xlabel('% Contrast')
% set(gca, 'XSCale', 'log')
title('CRF by FOV')

disp(['tetO:  n = ' num2str(size(dat,2)) ' by FOV'])
disp(['Ai203: n = ' num2str(size(dat2,2)) ' by FOV'])

%% save above data
% by FOV
tre_contrast = dat;
ai203_contrast = dat2;
save('e:/ai203/data/contrast_data_exported.mat', 'tre_contrast', 'ai203_contrast', '-v7.3')

%% save single cells
tre = crf;
ai203 = crf2;
save('e:/ai203/data/contrast_data_cells.mat', 'tre', 'ai203', '-v7.3')

%%

[All, outVars] = getTuningCurve(All, opts, outVars);
[All, outVars] = calcOSI(All, outVars);
[All, outVars] = calcTuningCircular(All, outVars); % note: only works on tuned cells (ie. not for max of visID=1)