function [stimID, uniqueStims, uStimCount] = stimReader(ExpStruct,epoch)
sws=[];
sws(1) = ExpStruct.EpochEnterSweep{epoch};
try
    sws(2) = ExpStruct.EpochEnterSweep{epoch+1}-1;
catch
    sws(2) = ExpStruct.sweep_counter-1;
end
sws = sws(1):sws(2);

stimID = ExpStruct.stim_tag(sws);

uniqueStims = unique(stimID);

uStimCount=[];
for i = 1:numel(uniqueStims)
    uStimCount(i)=sum(stimID==uniqueStims(i));
end
