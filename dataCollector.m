%% general info
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
disp('got info')

%% experiment epoch
exp.zdfData = zdfData;
exp.allData = allData;
exp.runVal = runVector;
exp.lowMotionTrials = lowMotionTrials;
exp.stimID = stimID;
exp.visID = visID;
exp.holoTargets = HoloTargets;
exp.rois = roisTargets;
exp.allCoM = allCoM;
exp.allDepth = allDepth;
exp.stimCoM = stimCoM;
exp.stimDepth = stimDepth;
exp.targetedCells = targettedCells;
exp.stimParams = stimParam;
exp.stimParams.Hz = holoRequests.holoStimParams.hzList;
exp.stimParams.numCells = holoRequests.holoStimParams.cellsPerHolo;
exp.stimParams.powers = holoRequests.holoStimParams.powerList;
disp('got exp')

%% orientation/vis epoch
vis.zdfData = zdfData;
vis.allData = allData;
vis.runVal = runVector;
vis.lowMotionTrials = lowMotionTrials;
vis.visID = visID;
disp('got vis')

%% run to save

out.info = info;
out.exp = exp;
out.vis = vis;

save([basePath info.date '_' info.mouse '_outfile'], 'out')
save(['Z:\willh\outputdata\' info.date '_' info.mouse 'outfile'], 'out')

disp('All data saved!')