
switch option
    case 'exp'
        if ~exist('ExpStruct')
            disp('Need to reload ExpStruct...')
            in = load(physfile,'ExpStruct');
            ExpStruct=in.ExpStruct;
        end
        
%         if length(date)>6
%             info.date = date(3:end);
%         else
            info.date = date;
%         end
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
        exp.dfData = dfData;
        
        exp.NPcoef = allNPC;
        exp.allNP = allNP;
        exp.allDataNoNP = allDataNoNP;
              
        
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
    case 'stim'
        %% stimTest
       
        stm.DAQepoch = DAQepoch;
        stm.zdfData = zdfData;
        stm.allData = allData;
        stm.dfData = dfData;
        
        stm.NPcoef = allNPC;
        stm.allNP = allNP;
        stm.allDataNoNP = allDataNoNP;
        
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
        
        disp('got stim')
        
    case 'vis'
        %% orientation/vis epoch
        vis.desc = 'Ori';
        vis.gratingSize = 50;
        vis.zdfData = zdfData;
        vis.allData = allData;
        vis.dfData = dfData;
        
        vis.NPcoef = allNPC;
        vis.allNP = allNP;
        vis.allDataNoNP = allDataNoNP;
        
        vis.runVal = runVector;
        vis.lowMotionTrials = lowMotionTrials;
        vis.visID = visID;
        vis.visStart = visStart;
        vis.visStop = visStop;
        vis.DAQepoch = DAQepoch;
        disp('got vis')
    case 'vis2'
        %% if two vis sets
        vis2.desc = 'GMN';
        vis2.zdfData = zdfData;
        vis2.allData = allData;
        vis2.dfData = dfData;
        
        vis2.NPcoef = allNPC;
        vis2.allNP = allNP;
        vis2.allDataNoNP = allDataNoNP;
        
        vis2.runVal = runVector;
        vis2.lowMotionTrials = lowMotionTrials;
        vis2.visID = visID;
        vis2.visStart = visStart;
        vis2.visStop = visStop;
        vis2.DAQepoch = DAQepoch;
        disp('got vis2')
    case 'vis3'
        %% if two vis sets
        vis3.desc = 'other';
        vis3.zdfData = zdfData;
        vis3.allData = allData;
        vis3.dfData = dfData;
        
        vis3.NPcoef = allNPC;
        vis3.allNP = allNP;
        vis3.allDataNoNP = allDataNoNP;
        
        vis3.runVal = runVector;
        vis3.lowMotionTrials = lowMotionTrials;
        vis3.visID = visID;
        vis3.visStart = visStart;
        vis3.visStop = visStop;
        vis3.DAQepoch = DAQepoch;
        disp('got vis3')
    case 'mani'
        %% save exp to Manifold
        mani = exp;
        try
            mani.mani = ExpStruct.Mani;
        catch
            try tempMani    = ExpStruct.ManifoldWrite; catch end
            tempMani.estSpikes      = ExpStruct.estSpikes;
            try
                tempMani.CellsWritten   = ExpStruct.CellsWritten;
                tempMani.rThreshold     = ExpStruct.rThreshold;
                tempMani.FitCells       = ExpStruct.FitCells;
            catch
                disp('************Do Not have all variables for Mani*******************')
            end
            mani.mani = tempMani;
        end
        try
            mani.CellIDs = unique(cat(2,exp.holoRequest.rois{:}));
        catch
            mani.CellIDs = unique(cat(1,exp.holoRequest.rois{:}));
        end
        
        try
            mani.estSpikes = ExpStruct.estSpikes;
        catch
        end
       case 'mani2'
        %% save exp to Manifold
        mani2 = exp;
        try
            mani2.mani = ExpStruct.Mani;
        catch
            try tempMani    = ExpStruct.ManifoldWrite; catch end
            tempMani.estSpikes      = ExpStruct.estSpikes;
            try
                tempMani.CellsWritten   = ExpStruct.CellsWritten;
                tempMani.rThreshold     = ExpStruct.rThreshold;
                tempMani.FitCells       = ExpStruct.FitCells;
            catch
                disp('************Do Not have all variables for Mani*******************')
            end
            mani2.mani = tempMani;
        end
        try
            mani2.CellIDs = unique(cat(2,exp.holoRequest.rois{:}));
        catch
            mani2.CellIDs = unique(cat(1,exp.holoRequest.rois{:}));
        end
        
        try
            mani2.estSpikes = ExpStruct.estSpikes;
        catch
        end     
        
    case 'spk'
        %% save exp to spk: for Spike Curve Support
        spk=exp;
        spk.dfData = dfData;
        try
            spk.Mani = ExpStruct.Mani;
        catch
            tempMani=[];
            try tempMani    = ExpStruct.spikeTest;catch;end
           
            try tempMani.estSpikes      = ExpStruct.estSpikes;catch;end
            
            try tempMani.rThreshold     = ExpStruct.rThreshold;catch;end
            try tempMani.FitCells       = ExpStruct.FitCells;catch;end
            try tempMani.pBalLookRight  = ExpStruct.pBalLookRight;catch;end
            try tempMani.roiWeightsPBal = ExpStruct.roiWeightsPBal;catch;end
            try tempMani.pBalSatValue   = ExpStruct.pBalSatValue;catch;end
            
                disp('****************Maybe Do not have all variables for spk****************')
            
            spk.Mani = tempMani;
        end
        disp('got spk')
        
    case 'exp2'
        %% experiment epoch
        exp2.DAQepoch = DAQepoch;
        exp2.zdfData = zdfData;
        exp2.allData = allData;
        exp2.dfData = dfData;
        
        exp2.NPcoef = allNPC;
        exp2.allNP = allNP;
        exp2.allDataNoNP = allDataNoNP;
        
        exp2.runVal = runVector;
        exp2.lowMotionTrials = lowMotionTrials;
        exp2.stimID = stimID;
        try
            exp2.visID = visID;
            exp2.visStart = visStart;
            exp2.visStop = visStop;
        catch;
            disp('No VisID')
        end
        exp2.holoTargets = HoloTargets;
        exp2.rois = roisTargets;
        exp2.allCoM = allCoM;
        exp2.allDepth = allDepth;
        exp2.stimCoM = stimCoM;
        exp2.stimDepth = stimDepth;
        exp2.targetedCells = targettedCells;
        exp2.uniqueStims = uniqueStims;
        exp2.outputsInfo = outputPatternTranslator(ExpStruct,uniqueStims);
        
        tempOutputOrder = exp2.outputsInfo.OutputOrder;
        tempOutputOrder(tempOutputOrder==0)=[];
        try exp2.output_names = ExpStruct.output_names; catch; end
        
        exp2.stimParams = stimParam;
        try
            exp2.stimParams.Hz = holoRequests.holoStimParams.hzList(tempOutputOrder);
            exp2.stimParams.numCells = holoRequests.holoStimParams.cellsPerHolo(tempOutputOrder);
            exp2.stimParams.powers = holoRequests.holoStimParams.powerList(tempOutputOrder);
        catch
            disp('no holoStimParams')
        end
        
        exp2.holoRequest = ExpStruct.Holo.holoRequests{...
            ExpStruct.Holo.Sweeps_holoRequestNumber(ExpStruct.EpochEnterSweep{DAQepoch})};
        
        exp2.Tarray = Tarray; %Motion Correct trace;
        exp2.dfData = dfData; %non zscored data;
        
        exp2.offsets = offsets; %sometimes different epochs calc subtly different offsets, recorded here
        
        disp('got exp')
    case 'info'
        %% general info
        
        if ~exist('ExpStruct')
            disp('Need to reload ExpStruct...')
            in = load(physfile,'ExpStruct');
            ExpStruct=in.ExpStruct;
        end
        
%         if length(date)>6
%             info.date = date(3:end);
%         else
            info.date = date;
%         end
        info.mouse = mouse;
        info.epochText1 = ExpStruct.EpochText1;
        info.epochText2 = ExpStruct.EpochText2;
        info.params = ExpStruct.Expt_Params;
        info.path = path;
        info.outputNames = ExpStruct.output_names;
        info.FR=FR;
        info.offsets = offsets; %motion correction and clipping offsets.
        
        disp('got info')
end
%%

try out.info = info; catch; end
try out.exp = exp; catch; end
try out.exp2 = exp2; catch; end;
try out.vis = vis; catch; end

try out.vis2=vis2; catch;end
try out.vis3=vis3; catch;end

try out.spk = spk; catch; end
try out.mani = mani; catch; end

try out.mani2 = mani2; catch; end
try out.mani0 = mani0; catch; end
try out.stm = stm; catch; end

out