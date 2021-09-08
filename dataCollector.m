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

%% experiment epoch
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
try exp.output_names = ExpStruct.output_names; catch; end

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

exp.offsets = offsets; %sometimes different epochs calc subtly different offsets, recorded here

disp('got exp')

%% save exp to spk: for Spike Curve Support
spk=exp;
spk.dfData = dfData;
try
spk.Mani = ExpStruct.Mani;
catch
    tempMani    = ExpStruct.spikeTest;
    tempMani.estSpikes      = ExpStruct.estSpikes;
    tempMani.rThreshold     = ExpStruct.rThreshold;
    tempMani.FitCells       = ExpStruct.FitCells;
    tempMani.pBalLookRight  = ExpStruct.pBalLookRight;
    tempMani.roiWeightsPBal = ExpStruct.roiWeightsPBal;
    tempMani.pBalSatValue   = ExpStruct.pBalSatValue;
    spk.Mani = tempMani;
end
disp('got spk')

%% save exp to Manifold
mani = exp;
try
mani.mani = ExpStruct.Mani;
catch
    tempMani    = ExpStruct.ManifoldWrite;
    tempMani.CellsWritten   = ExpStruct.CellsWritten;
    tempMani.estSpikes      = ExpStruct.estSpikes;
    tempMani.rThreshold     = ExpStruct.rThreshold;
    tempMani.FitCells       = ExpStruct.FitCells;
    mani.mani = tempMani;
end
mani.CellIDs = unique(cat(2,exp.holoRequest.rois{:}));

%% orientation/vis epoch
vis.desc = 'Ori';
vis.gratingSize = 50;
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


%% stimTest
stm.DAQepoch = DAQepoch;
stm.zdfData = zdfData;
stm.allData = allData;
stm.stimID =stimID;

swpStart = ExpStruct.EpochEnterSweep{DAQepoch};
Hnum = ExpStruct.Holo.Sweeps_holoRequestNumber(swpStart);
HR =ExpStruct.Holo.holoRequests{Hnum-1};
HR2 =ExpStruct.Holo.holoRequests{Hnum};

stm.holoRequest = HR;
stm.holoRequest2 = HR2; %sometimes its putting the wrong HR in save both.

stm.runVal = runVector;
stm.lowMotionTrials = lowMotionTrials;

stm.holoTargets = HoloTargets;
stm.rois = roisTargets;
stm.allCoM = allCoM;
stm.allDepth = allDepth;
stm.stimCoM = stimCoM;
stm.stimDepth = stimDepth;
stm.targetedCells = targettedCells;
stm.uniqueStims = uniqueStims;
stm.outputsInfo = outputPatternTranslator(ExpStruct,uniqueStims);

tempOutputOrder = stm.outputsInfo.OutputOrder;
tempOutputOrder(tempOutputOrder==0)=[];
% exp.output_names = ExpStruct.output_names;

stm.stimParams = stimParam;
try
    stm.stimParams.Hz = holoRequests.holoStimParams.hzList(tempOutputOrder);
    stm.stimParams.numCells = holoRequests.holoStimParams.cellsPerHolo(tempOutputOrder);
    stm.stimParams.powers = holoRequests.holoStimParams.powerList(tempOutputOrder);
catch
    try
        stm.stimParams.Hz = holoRequests.holoStimParams.hzList(1);
        stm.stimParams.numCells = holoRequests.holoStimParams.cellsPerHolo(1);
        stm.stimParams.powers = holoRequests.holoStimParams.powerList(1);
        disp('holostimparams only one value used')
    catch
        disp('no holoStimParams')
    end
end

stm.Tarray = Tarray; %Motion Correct trace;
stm.dfData = dfData; %non zscored data;

stm.offsets = offsets; %sometimes different epochs calc subtly different offsets, recorded here


%% run to save

out.info = info;
out.exp = exp;
try %if merging two or more experiments put them in here
    out.exp1 = exp1;
    out.exp2 = exp2;
catch;
end;
try
    out.vis = vis;
catch
    disp('No Vis Data')
end

try out.vis2=vis2; catch;end
try out.spk = spk; catch; end
try out.mani = mani; catch; end
try out.mani2 = mani2; catch; end
try out.mani0 = mani0; catch; end
try out.stm = stm; catch; end

out

save([basePath info.date '_' info.mouse '_outfile'], 'out','-v7.3')
% save(['Z:\willh\outputdata\' info.date '_' info.mouse 'outfile'], 'out')
% save(['U:\ioldenburg\outputdata1\' info.date '_' info.mouse '_outfile'], 'out')
% save(['E:\Contrast Modulated Ensembles\'

disp('All data saved!')

%% Merge exps. simple experiments that can be merged should be done so here
%be aware that some identical experiments, can have different numbers of
%frames

ex1 = exp1;exp_e5;% exp1;
ex2 = exp_e9;% exp2;

ex=ex1;

frameLen = min(size(ex1.zdfData,2),size(ex2.zdfData,2));

ex.zdfData = cat(3,ex1.zdfData(:,1:frameLen,:),ex2.zdfData(:,1:frameLen,:));
ex.allData = cat(3,ex1.allData(:,1:frameLen,:),ex2.allData(:,1:frameLen,:));
ex.runVal = cat(1,ex1.runVal,ex2.runVal);
ex.lowMotionTrials = cat(2,ex1.lowMotionTrials,ex2.lowMotionTrials);
ex.stimID = cat(2,ex1.stimID,ex2.stimID);
ex.visID = cat(2,ex1.visID,ex2.visID);
ex.DAQepoch = [ex1.DAQepoch ex2.DAQepoch];
ex.dfData = cat(3,ex1.dfData(:,1:frameLen,:),ex2.dfData(:,1:frameLen,:));

ex.uniqueStims = unique([ex1.uniqueStims ex2.uniqueStims]);
ex.holoTargets = [ex1.holoTargets ex2.holoTargets];
ex.rois = [ex1.rois ex2.rois];

numRois = numel(ex1.rois); %warning assumes first seq of ex2 is 0, ie no stim.
Seq = [ex1.stimParams.Seq ex2.stimParams.Seq(2:end)+numRois];
numPulse = [ex1.stimParams.numPulse ex2.stimParams.numPulse(2:end)];
roi = [ex1.stimParams.roi cellfun(@(x) x+numRois, ex2.stimParams.roi(2:end),'uniformoutput',0)];
Hz = [ex1.stimParams.Hz ex2.stimParams.Hz];
numCells = [ex1.stimParams.numCells ex2.stimParams.numCells];
powers = [ex1.stimParams.powers ex2.stimParams.powers];

stimParams.Seq = Seq;
stimParams.numPulse = numPulse;
stimParams.roi = roi;
stimParams.Hz = Hz;
stimParams.numCells = numCells;
stimParams.powers = powers;

ex.stimParams = stimParams;


oi.OutputStims = [ex1.outputsInfo.OutputStims ex1.outputsInfo.OutputStims(2:end)];
oi.OutputNames = [ex1.outputsInfo.OutputNames ex1.outputsInfo.OutputNames(2:end)];
oi.OutputOrder = [ex1.outputsInfo.OutputOrder ex1.outputsInfo.OutputOrder(2:end)];
oi.OutputPatterns = [ex1.outputsInfo.OutputPatterns ex1.outputsInfo.OutputPatterns(2:end)];
ex.outputsInfo = oi;

for i =1:numel(ex1.Tarray)
    ex.Tarray{i} = cat(2,ex1.Tarray{i}(1:frameLen,:,:),ex2.Tarray{i}(1:frameLen,:,:));
end

exp = ex;
