oriLoadList;
loadPath = 'T:\Outfiles';

addpath(genpath(loadPath))

numExps = numel(loadList);

%% Load data

numExps = numel(loadList);
disp(['There are ' num2str(numExps) ' Exps in this LoadList'])
if numExps ~= 0
    clear All
    if ~iscell(loadList)
        numExps=1;
        temp = loadList;
        clear loadList;
        loadList{1} = temp;
    end
    for ind = 1:numExps
        pTime =tic;
        fprintf(['Loading Experiment ' num2str(ind) '...']);
        All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
        fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
    end
else
    disp('Did you press this by accident?')
end


%% add neuropil and neuropil coefficient to out
savePath = [loadPath '\new'];

willexps = [11 14 15 16 19 20 21 22 23 24 25 29];

for ind=willexps; %30:numExps;
    dataPathBackup=[];
    if ~isfield(All(ind).out.exp,'allNP');
        pTime =tic;
        fprintf(['Loading Experiment ' num2str(ind) '...']);
        
        in = load(fullfile(loadPath,loadList{ind}),'out');
        out = in.out;
        
        %get NP coefficient
        allDepth = out.exp.allDepth;
        numCells = size(out.exp.allData,1);
        NPcoef = zeros([1 numCells]);
        for cellID = 1:numCells
            depth = allDepth(cellID);
            
            [dat, All(ind)] = pullDat(All(ind),depth);
            
            cellList = 1:numel(allDepth);
            %figure out the ID in as dat sees it
            IDinDat = numel(find(allDepth==depth & cellList'<=cellID));
            
            iscelldepth = [dat.stat.iscell];
            iscelldepth = find(iscelldepth);
            IDinDat = iscelldepth(IDinDat);
            
            NPcoef(cellID) = dat.stat(IDinDat).neuropilCoefficient;
        end
        out.exp.NPcoef = NPcoef;
        disp('got NP coefs')
        
        %% NP Values
        disp('Working on NP values')
        nDepthsTotal =   max(allDepth);
        numColors = 2;
        
        %get dat
        disp('getting dat')
        %         eval(['dat = All(ind).out.dat' num2str(d) ';'])
        dat = All(ind).out.dat1;
        
        %determine DAQepoch to work from the right section
        if isfield(out.exp,'DAQepoch')
            DAQepoch = out.exp.DAQepoch;
        else
            disp('EpochText 1')
            All(ind).out.info.epochText1
            disp('Epoch Text 2')
            All(ind).out.info.epochText2
            expNumFiles = size(All(ind).out.exp.allData,3);
            disp(['Expect ' num2str(expNumFiles) ' Files.']);
            disp([All(ind).out.info.mouse ' ' All(ind).out.info.date])
            e = input('What was the DAQepoch Num?');
            DAQepoch = e;
            out.exp.DAQepoch = e;
        end
        
        rootpath = dat.ops.RootDir;
        dataPath = fullfile(rootpath,num2str(DAQepoch));
        
        if ~exist(dataPath) && isempty(dataPathBackup) && ( ~isfield(out.info,'dataPath2') || isempty(out.info.dataPath2))
            disp('Path Not Found')
            disp(dataPath)
            dataPath = uigetdir('Y:\frankenshare\ian','Select Data Path');
%             dataPath = uigetdir('Set_best_guess_path','Select Data Path');
            dataPathBackup = dataPath;
            out.info.dataPath2 = dataPath; 
            
        elseif isfield(out.info,'dataPath2') && ~isempty(out.info.dataPath2)
            disp('using out.info.dataPath2')
            dataPath = out.info.dataPath2;
        elseif ~isempty(dataPathBackup)
            disp('pulling path from backup')
            dataPath = dataPathBackup;
        end
        
        fprintf('Identifying Files... \n');
        k=dir(dataPath);k(1:2)=[];
        localFiles=[];
        i=1;
        for n=1:(numel(k))
            fn=k(n).name;
            if regexp(fn,regexptranslate('wildcard',['*.tif']))
                localFiles{i}=fullfile(dataPath,  k(n).name);
                i=i+1;
            end;
        end;
        fprintf([ num2str(numel(localFiles)) ' Files Detected... \n']);
        
        %%Determine number of frames per tiff
        fcTime = tic;
        try
            temp = load(fullfile(dataPath,'FC.mat'));
        catch
            temp=[];
        end
        
        if isfield(temp,'frameCount')
            disp('Saved Frame Count Detected, skipping load')
            frameCount=temp.frameCount;
        else
            disp('No Frame Count Detected Loading Files...')
            
            frameCount=[];
            for i = 1:numel(localFiles)
                fprintf(['Frame ' num2str(i) ' of ' num2str(numel(localFiles)) '. ']);
                tic;
                try
                    x = ScanImageTiffReader(localFiles{i}).descriptions();
                    frameCount(i) = size(x,1)/nDepthsTotal/numColors; %assumes two color
                catch
                    fprintf(['Used iminfo... ']);
                    im = imfinfo(localFiles{i});
                    frameCount(i) = numel(im)/numColors/nDepthsTotal;
                end
                fprintf([num2str(toc) 's\n']);
            end
            
            disp('Saving frameCount')
            save(fullfile(dataPath,'FC.mat'),'frameCount')
            
            
            
            toc(fcTime)
        end
        
        %% extract the NP data 
        
        minNumFrames=min(frameCount);
        numFrames =floor(minNumFrames);
        
        s2pEpoch =find(dat.ops.expts==DAQepoch);
        
         allDataNoNP=[];  allNP=[];
         clear FdataNew FdataNoNPNew FNPNew
        for d=1:max(allDepth);
            eval(['dat = All(ind).out.dat' num2str(d) ';'])
            
            cellID = [dat.stat(:).iscell];
            
            temp = dat.Fcell{s2pEpoch};
            Fdata = temp(find(cellID'==1),:);
            Fnp = dat.FcellNeu{s2pEpoch}(find(cellID'==1),:);
            npc = [dat.stat(find(cellID'==1)).neuropilCoefficient];
            FdataNoNP = Fdata;
            %             Fdata= Fdata-Fnp.*npc';
            
            sz=size(Fdata);
            
            FdataNew=[];FdataNoNPNew=[];
            for i=1:numel(frameCount)
                FdataNew{i} = Fdata(:,sum(frameCount(1:i-1))+1 : min(sum(frameCount(1:i)),sz(2)) );
                FdataNoNPNew{i} = FdataNoNP(:,sum(frameCount(1:i-1))+1 : min(sum(frameCount(1:i)),sz(2)) );
                FNPNew{i} = Fnp(:,sum(frameCount(1:i-1))+1 : min(sum(frameCount(1:i)),sz(2)) );
            end
            
            temp=[];
            temp=cellfun(@(x) x(:,1:numFrames),FdataNoNPNew,'UniformOutput',false);
            FdataNoNP = cell2mat(temp);
            FdataNoNP = reshape(FdataNoNP,sz(1),numFrames,numel(frameCount));
            
            temp=[];
            temp=cellfun(@(x) x(:,1:numFrames),FNPNew,'UniformOutput',false);
            Fnp = cell2mat(temp);
            Fnp = reshape(Fnp,sz(1),numFrames,numel(frameCount));
            
            allDataNoNP = cat(1,allDataNoNP,FdataNoNP);
            allNP = cat(1,allNP,Fnp);            
        end
        out.exp.allNP = allNP;
        out.exp.allDataNoNP = allDataNoNP;
        
        
        
        %%
        save(fullfile(savePath,loadList{ind}),'out','-v7.3')
        fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
    end
end
disp('done')