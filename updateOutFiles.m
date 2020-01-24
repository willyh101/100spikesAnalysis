[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

%% load
numExps = numel(loadList);

clear All
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% ID lack of offsetts

for ind = 1:numExps
    if ~isfield(All(ind).out.info,'offsets')
         pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    try
        out = All(ind).out;
        pth = fileparts(All(ind).out.info.path);
        x = dir(fullfile(pth,'Analysis*'));
        
        inDat = load(fullfile(pth,x(1).name),'offsets','FR');
        out.info.offsets =  inDat.offsets;
        fprintf(['\nOffsets were: ' num2str(inDat.offsets) '\n']);
        out.info.FR  = inDat.FR;
        
        save(fullfile(loadPath,loadList{ind}),'out')
    catch
        disp(['Unable to edit Experiment ' num2str(ind)])
    end
        
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

    end
    
end
disp('Done')


%%
%% ID lack of visStart

for ind = 1:numExps
    if ~isfield(All(ind).out.exp,'visStart')
         pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    try
        out = All(ind).out;
        pth = fileparts(All(ind).out.info.path);
        x = dir(fullfile(pth,'Analysis*'));
        
        inDat = load(fullfile(pth,x(1).name),'visStart','visStop');
       out.exp.visStart = inDat.visStart;
       out.exp.visStop = inDat.visStop;
        
        fprintf(['\nVisStart was: ' num2str(inDat.visStart) '\n']);
        
        save(fullfile(loadPath,loadList{ind}),'out')
    catch
        disp(['Unable to edit Experiment ' num2str(ind)])
    end
        
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

    end
    
end
disp('Done')

%% ID and fix lack of outputPaterns

for ind = 1:numExps
    if ~isfield(All(ind).out.exp,'outputsInfo')
         pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    try
        out = All(ind).out;
        pth = fileparts(All(ind).out.info.path);
        x = dir(fullfile(pth,'Analysis*'));
        
        inDat = load(fullfile(pth,x(1).name),'physfile');
        inDat2 = load(inDat.physfile,'ExpStruct');
        
        uniqueStims =unique(All(1).out.exp.stimID);
        out.exp.outputsInfo = outputPatternTranslator(inDat2.ExpStruct,uniqueStims);
       
             
        save(fullfile(loadPath,loadList{ind}),'out')
    catch
        disp(['Unable to edit Experiment ' num2str(ind)])
    end
        
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

    end
    
end
disp('Done')
%% add vis stim IDs to out files
ptime=tic;
clear invar
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    invar(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

for ind = 1:numExps
    out = invar(ind).out;
    vises = All(ind).out.exp.visCond;
    out.exp.visCond = vises;
    save(fullfile(savePath,loadList{ind}),'out')
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])   
end
disp('Done')







