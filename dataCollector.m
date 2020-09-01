%% general info

if ~exist('ExpStruct')
    disp('Need to reload ExpStruct...')
in = load(physfile,'ExpStruct');
ExpStruct=in.ExpStruct;
end

if length(date)>6
    info.date = date(3:end);
else
    info.date = date;
end
info.mouse = mouse;
info.epochText1 = ExpStruct.EpochText1;
info.epochText2 = ExpStruct.EpochText2;
info.params = ExpStruct.Expt_Params;
info.path = path;
info.outputNames = ExpStruct.output_names;
info.FR=FR;
info.offsets = offsets; %motion correction and clipping offsets.

disp('got info')

%%experiment epoch
exp.DAQepoch = DAQepoch;
exp.zdfData = zdfData;
exp.allData = allData;
exp.runVal = runVector;
exp.lowMotionTrials = lowMotionTrials;
exp.stimID = stimID;
try
exp.visID = visID;
exp.visStart = visStart;
exp.visStop = visStop;
catch;
disp('No VisID')
end
exp.holoTargets = HoloTargets;
exp.rois = roisTargets;
exp.allCoM = allCoM;
exp.allDepth = allDepth;
exp.stimCoM = stimCoM;
exp.stimDepth = stimDepth;
exp.targetedCells = targettedCells;
exp.uniqueStims = uniqueStims;
exp.outputsInfo = outputPatternTranslator(ExpStruct,uniqueStims);

tempOutputOrder = exp.outputsInfo.OutputOrder;
tempOutputOrder(tempOutputOrder==0)=[];
% exp.output_names = ExpStruct.output_names;

exp.stimParams = stimParam;
try
exp.stimParams.Hz = holoRequests.holoStimParams.hzList(tempOutputOrder);
exp.stimParams.numCells = holoRequests.holoStimParams.cellsPerHolo(tempOutputOrder);
exp.stimParams.powers = holoRequests.holoStimParams.powerList(tempOutputOrder);
catch
    disp('no holoStimParams')
end

exp.holoRequest = ExpStruct.Holo.holoRequests{...
    ExpStruct.Holo.Sweeps_holoRequestNumber(ExpStruct.EpochEnterSweep{DAQepoch})};

exp.Tarray = Tarray; %Motion Correct trace;
exp.dfData = dfData; %non zscored data; 

disp('got exp')

%% save exp to spk: for Spike Curve Support
spk=exp;
spk.dfData = dfData;
spk.Mani = ExpStruct.Mani;

%% save exp to Manifold
mani = exp;
mani.mani = ExpStruct.Mani; 

%% orientation/vis epoch
vis.desc = 'Ori';
vis.zdfData = zdfData;
vis.allData = allData;
vis.runVal = runVector;
vis.lowMotionTrials = lowMotionTrials;
vis.visID = visID;
vis.visStart = visStart;
vis.visStop = visStop;
vis.DAQepoch = DAQepoch;
disp('got vis')

%% if two vis sets
vis2.desc = 'GMN';
vis2.zdfData = zdfData;
vis2.allData = allData;
vis2.runVal = runVector;
vis2.lowMotionTrials = lowMotionTrials;
vis2.visID = visID;
vis2.visStart = visStart;
vis2.visStop = visStop;
vis2.DAQepoch = DAQepoch;
disp('got vis2')
%% run to save

out.info = info;
out.exp = exp;
try %if merging two or more experiments put them in here
out.exp1 = exp1;
out.exp2 = exp2;
catch; end;
try
    out.vis = vis;
catch
    disp('No Vis Data')
end

save([basePath info.date '_' info.mouse '_outfile'], 'out')
% save(['Z:\willh\outputdata\' info.date '_' info.mouse 'outfile'], 'out')
% save(['U:\ioldenburg\outputdata1\' info.date '_' info.mouse '_outfile'], 'out')
% save(['E:\Contrast Modulated Ensembles\'

disp('All data saved!')

%% Merge exps. simple experiments that can be merged should be done so here
%be aware that some identical experiments, can have different numbers of
%frames

ex1 = exp1;
ex2 = exp2;

ex=ex1;

frameLen = min(size(ex1.zdfData,2),size(ex2.zdfData,2));

ex.zdfData = cat(3,ex1.zdfData(:,1:frameLen,:),ex2.zdfData(:,1:frameLen,:));
ex.allData = cat(3,ex1.allData(:,1:frameLen,:),ex2.allData(:,1:frameLen,:));
ex.runVal = cat(1,ex1.runVal,ex2.runVal);
ex.lowMotionTrials = cat(2,ex1.lowMotionTrials,ex2.lowMotionTrials);
ex.stimID = cat(2,ex1.stimID,ex2.stimID);
ex.visID = cat(2,ex1.visID,ex2.visID);
ex.DAQepoch = [ex1.DAQepoch ex2.DAQepoch];

exp = ex;
