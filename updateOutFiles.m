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
