ai203LoadList;
treChroMELoadList
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

%% Check if stm is present and if not add

for ind = 1:numExps
    if isfield(All(ind).out, 'stm')
        disp(['ind: ' num2str(ind) ' stim found.'])
    else
        
        %
        out = All(ind).out;
        info = All(ind).out.info;
        
        disp(['Processing: ' info.date '_' info.mouse]);
        
        fullPth = info.path;
        [pth, epochNums] = fileparts(info.path);
        date = info.date;
        
        inputParams=[];
        inputParams.fullPth = fullPth;
        inputParams.date = date;
        
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
        
        inputParams.physfile = fullfile(pth,physfileName);
        if ~isempty(dir(inputParams.physfile))
            
            baseName = [info.mouse '_' info.date];
            loadListPlanes = {['F_' baseName '_plane1_proc'] ['F_' baseName '_plane2_proc'] ['F_' baseName '_plane3_proc']};% ['F_' baseName '_plane4_proc']};
            inputParams.loadList = loadListPlanes;
            
            %select epoch to export
            info.epochText1(:)
            y = input('Select epoch to use for stim, 0  launch Acq, or anything else to abort ','s');
            y = str2num(y);
            
            if  y==0
                load(inputParams.physfile);
                disp('Loaded physfile');
                return
            elseif ~isnumeric(y) || isempty(y) || y<0
                disp('Aborting')
                return
            else
                disp(['you selected Epoch ' num2str(y) '. "' info.epochText1{y} '"'])
                inputParams.DAQepoch = y;
                disp(['S2P Epochs: ' epochNums]);
                s2pEpoch = find(cellfun(@(x) str2num(x)==y,split(epochNums,'_')));
                disp(['Detected s2p epoch is ' num2str(s2pEpoch)]);
                inputParams.s2pEpoch = s2pEpoch;
            end
            
%             key = input('Ready to Procede? (y/n)','s');
%             if ~strcmp(key,'y')
%                 disp('Aborted')
%                 return
%             else
                disp('Launching...')
                stm = ExportStm(inputParams); %Main Function!
                out.stm = stm;
                
                disp('Saving...')
                save(fullfile(savePath,loadList{ind}),'out')
                disp('done')
%             end
        else
            disp('Could not load file, maybe wrong directory')
            disp(inputParams.physfile)
        end
        
    end
end
disp('All Done')