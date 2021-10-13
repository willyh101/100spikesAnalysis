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
        
        uniqueStims =unique(All(ind).out.exp.stimID);
        if any(uniqueStims==0)
            disp('Unholly abomination of a stimID==0. Be afraid')
            uniqueStims(uniqueStims==0)=[];
        end
        
        out.exp.outputsInfo = outputPatternTranslator(inDat2.ExpStruct,uniqueStims);
       
             
        save(fullfile(loadPath,loadList{ind}),'out')
    catch
        disp(['Unable to edit Experiment ' num2str(ind)])
    end
        
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

    end
    
end
disp('Done')
%%



%% load will's out files

for ind = [15 16]
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    try
        out = All(ind).out;
        pth = fileparts(All(ind).out.info.path);
        inDat = load(fullfile(pth, 'physfile.mat'), 'ExpStruct');
        %inDat2 = load(inDat.physfile,'ExpStruct');
        uniqueStims =unique(All(ind).out.exp.stimID);
        out.exp.outputsInfo = outputPatternTranslator(inDat.ExpStruct,uniqueStims);
        save(fullfile(savePath,loadList{ind}),'out')
    catch
        disp(['Unable to edit Experiment ' num2str(ind)])
    end
        
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

end
    

disp('Done')
%% add vis stim IDs to out files

%% ID vis stim conditions
% load in visConds.mat
clear visCond
for i=1:numel(All)
    All(i).out.exp.visKey = vis_out{i};
    visCond = [];
    if isempty(All(i).out.exp.visKey)
        continue
    end

    for k=1:numel(All(i).out.exp.visID)
        if isempty(All(i).out.exp.visKey)
            visCond(:,k)=nan(2,1);
            continue
        else
            visCond(:,k) = vis_out{i}(:,All(i).out.exp.visID(k));
        end
    end

    All(i).out.exp.visCond = visCond;
end
%%
% test = All(10).out.exp.visID(i);
% 
% arr = [0:0.2:1];
% for i=1:numel(All(10).out.exp.visID)
%     vis=All(10).out.exp.visID(i);
%     idx = find(arr==vis);
%     tempVis(i) = idx;
% end
% All(10).out.exp.visID = tempVis;





%%
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




%% Compute and add .dfData to outfiles
savePath = [loadPath '\new'];

for ind=1:numExps;
    if ~isfield(All(ind).out.exp,'dfData');
        pTime =tic;
        fprintf(['Loading Experiment ' num2str(ind) '...']);
        
        in = load(fullfile(loadPath,loadList{ind}),'out');
        out = in.out;
        [dfData, ~] =  computeDFFwithMovingBaseline(out.exp.allData);
        out.exp.dfData = dfData;
        
        
        save(fullfile(savePath,loadList{ind}),'out')
        fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
    end
end
disp('done')

%%

% ind = 55;
%  pTime =tic;
%         fprintf(['Loading Experiment ' num2str(ind) '...']);
%         
%         in = load(fullfile(loadPath,loadList{ind}),'out');
%         out = in.out;
%         out.exp = out.exp1;
%         
% %         [dfData, ~] =  computeDFFwithMovingBaseline(out.exp.allData);
% %         out.exp.dfData = dfData;
%         
%         
%         save(fullfile(savePath,loadList{ind}),'out')
%         fprintf([' Took ' num2str(toc(pTime)) 's.\n'])