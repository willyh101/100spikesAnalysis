function [visTable] = makeVisTable(out)

names = {'zdf', 'visID', 'lowMotion', 'lowRun', 'highRun'};
dataIn = out.vis;

rdata = dataIn.rdata';
vis = dataIn.visID';
lowMotion = dataIn.lowMotionTrials';
lowRun = dataIn.lowRunTrials';
highRun = ~lowRun;

visTable = table(rdata, vis, logical(lowMotion), logical(lowRun),...
    logical(highRun), 'VariableNames', names);



