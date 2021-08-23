ai203LoadList;
% testLoadList
loadPath = 'T:\Outfiles';
savePath = 'T:\Outfiles\temp';

%% load
numExps = numel(loadList);

clear All
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end
disp('Loaded')

%%
for ind = 1:numExps
    if isfield(All(ind).out, 'stm')
        disp(['ind: ' num2str(ind) ' stim found.'])
        
        out = All(ind).out;
        
        HR = All(ind).out.stm.holoRequest;
        stm = All(ind).out.stm;
        date = All(ind).out.info.date;
        
        if ~isfield(HR,'bigListOfFirstStimTimes')
            info = All(ind).out.info;
            [pth, epochNums] = fileparts(info.path);
            physfileName=[];
            x=dir(pth);
            for i =1:numel(x)
                if strcmp(x(i).name,[date '_A.mat'])...
                        || strcmp(x(i).name,[date '_B.mat'])...
                        || strcmp(x(i).name,[date '_C.mat']) ;
                    physfileName = x(i).name;
                end
            end
            if isempty(physfileName)
                disp('Error IDing phys file')
            end
            disp('Loading physfile');

            inputParams.physfile = fullfile(pth,physfileName);
            in = load(inputParams.physfile,'ExpStruct');
            ExpStruct = in.ExpStruct;
            
            DAQepoch = All(ind).out.stm.DAQepoch;
            swpStart = ExpStruct.EpochEnterSweep{DAQepoch};
            Hnum = ExpStruct.Holo.Sweeps_holoRequestNumber(swpStart);
            HR =ExpStruct.Holo.holoRequests{Hnum-1};
            
            stm.holoRequest = HR;
            
            out.stm =stm;
            
            disp('Saving Outfile')
            disp('Saving...')
            save(fullfile(savePath,loadList{ind}),'out')
            disp('done')
        end
    end
end


            
            