% outfile readme
% 
% out.info
% 
% .date           - date
% .mouse          - mouse name
% .epochText1     - description of each epoch
% .epochText2     - description of each epoch cont
% .params         - Experimental Notes includes age, imaging parameters, and more
% .path           - data path
% .outputNames    - holographic stim condition names. 
% .FR             - Frame Rate
% .offsetts       - registration offsets xy
% 
% out.exp
% .DAQepoch       - Epoch used 
% .zdfData        - Processed data zscored df in the order [cells frames trials]
% .dfData         - no zscore just df data
% .allData        - Raw Data
% .runVal         - run speed in cm/s 
% .lowMotionTrials - binary if low motion during stim period. 
% 
% .stimID         - hologram trial type. (the number is arbitrary usually indexes into out.info.outputNames but errors can occur)
% .uniqueStims    - Unique stimIDs. 
% .outputsInfo    - out of order control code
%     .OutputOrder    - if this is not monotonic you have an out of order error and uniquestims can be resorted using this field. formally indexes into stimParams
%     .OutputStims    - index of uniqueStims into outputStims points to the correct name in out.info.outputNames
%     .OutputName     - another place where the name is
%     .OutputPatterns - the actuall waveforms of command to the laser / vis computer / SI computer for that stim type
% .stimParams     - code for understanding stims. uniquestims position corresponds to 
%     .Seq        - the hologram identier. pulses sent to holocomputer to dictate which holgram used. 
%     .numPulses  - the number of laser pulses
%     .roi        - the roi identifier used in that hologram (see below)
%     .Hz         - rate (0 condition omited)
%     .numCells   - numCells in that hologram (0 omited)
%     .powers     - the nominal power per target used (note: this is transformed by the stimmability and referenced to 50mW so will be higher)
% 
% .visID          - the vis stimulus that trial. (0 is no feedback from vis computer, the rest of the code is in order set by vis computer)
% .visStart       - when vis started
% .visStop        - when vis ended
% 
% .rois           - cell array of the targets shot in a given hologram (target = makemask3d index).  index of cell array is set by out.exp.stimParam.roi
% .holoTargets    - as above but translated to matched cells as in .zdfData. (NaNs did not match)
% .targetedCells  - the match onlineData to matched Cells
% 
% .allCoM         - xy location of cells
% .allDepth       - z location of cells
% .stimCoM        - the xy location of targets
% .stimDepth      - the x location of targets
% 
% .holoRequest    - saved data sent to holostim computer
% .Tarray         - the full motion correction dataset
% 
% 
% out.vis
% 
% .desc           - type of vis experiment ori or contrast
% .zdfData
% .allData
% .runVal
% .lowMotionTrials
% .visID
% .visStart
% .visStop
% .DAQepoch       - all as above
% 
% out.Red 
% 
% to be continued
