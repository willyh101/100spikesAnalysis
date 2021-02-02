%data collector no daq

if length(date)>6
    info.date = date(3:end);
else
    info.date = date;
end
info.mouse = mouse;
info.path = path;
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

exp.holoTargets = HoloTargets;
exp.rois = roisTargets;
exp.allCoM = allCoM;
exp.allDepth = allDepth;
exp.stimCoM = stimCoM;
exp.stimDepth = stimDepth;
exp.targetedCells = targettedCells;
exp.uniqueStims = uniqueStims;

exp.Tarray = Tarray; %Motion Correct trace;
exp.dfData = dfData; %non zscored data; 

exp.holoRequest = holoRequest; 

exp.visID = visID;


%% %% orientation/vis epoch
vis.desc = 'Ori';
vis.zdfData = zdfData;
vis.allData = allData;
vis.runVal = runVector;
vis.lowMotionTrials = lowMotionTrials;
vis.visID = visID;
% vis.visStart = visStart;
% vis.visStop = visStop;
vis.DAQepoch = DAQepoch;
disp('got vis')