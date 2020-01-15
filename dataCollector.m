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
exp.zdfData = zdfData;
exp.allData = allData;
exp.runVal = runVector;
exp.lowMotionTrials = lowMotionTrials;
exp.stimID = stimID;
try
exp.visID = visID;
catch;end
exp.holoTargets = HoloTargets;
exp.rois = roisTargets;
exp.allCoM = allCoM;
exp.allDepth = allDepth;
exp.stimCoM = stimCoM;
exp.stimDepth = stimDepth;
exp.targetedCells = targettedCells;
exp.outputsInfo = outputPatternTranslator(ExpStruct,uniqueStims);

exp.stimParams = stimParam;
try
exp.stimParams.Hz = holoRequests.holoStimParams.hzList;
exp.stimParams.numCells = holoRequests.holoStimParams.cellsPerHolo;
exp.stimParams.powers = holoRequests.holoStimParams.powerList;
catch
    disp('no holoStimParams')
end
disp('got exp')

%% orientation/vis epoch
vis.zdfData = zdfData;
vis.allData = allData;
vis.runVal = runVector;
vis.lowMotionTrials = lowMotionTrials;
vis.visID = visID;
vis.visStart = visStart;
vis.visStop = visStop; 
disp('got vis')

%% run to save

out.info = info;
out.exp = exp;
try
    out.vis = vis;
catch
    disp('No Vis Data')
end

save([basePath info.date '_' info.mouse '_outfile'], 'out')
save(['Z:\willh\outputdata\' info.date '_' info.mouse 'outfile'], 'out')

disp('All data saved!')