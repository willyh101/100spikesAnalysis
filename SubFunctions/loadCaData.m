function   [localFiles, frameCount, MD, FR, minNumFrames, numFrames, numCells,...
    numFiles, allData, allSPData, allCoM, allDepth, allNP, allNPC,...
    roiMasks, T, BL, allDataCorrect, zData, zSPData, backupData,...
    DFFdata, backupSPData, ExpStruct, CellstToExclude, allDataNoNP,...
    stimID, uniqueStims, uStimCount, runVector, ArtifactSizeLeft,...
    ArtifactSizeRight, ROIinArtifact, excludeCells, framesUsed, Tused, zData2]...
    = loadCaData(path,loadList,DAQepoch,s2pEpoch,nDepthsTotal,ExpStruct,DIGITSWEEP,BaseLinePeriod,numColors)


UseNPCorrection=1;


load(fullfile(path,[loadList{1} '.mat']),'dat');
rootpath = dat.ops.RootDir;
dataPath = fullfile(rootpath,num2str(DAQepoch));

if ~exist(dataPath)
    dataPath = uigetdir('C:\Users\ian\Documents\DATA\DataStorage\','Select Data Path');
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
    fprintf(['Frame ' num2str(i) '. ']);
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
end

toc(fcTime)

localFiles(frameCount==1)=[];
frameCount(frameCount==1)=[];

if ~exist('im')
    im = imfinfo(localFiles{1});
end

MD = parseSI5Header(im(1).Software);
FR = MD.SI.hRoiManager.scanVolumeRate;
minNumFrames=min(frameCount);
% if ~exist('numFrames')
numFrames =floor(minNumFrames);
% end

numFiles = numel(frameCount);

%%True Loading Part
allData =[];
allDataNoNP=[];
allSPData=[];
allCoM=[];
allDepth =[];
tempCOM=[];
allNP = [];
allNPC = [];

skipSP = 0;

% Mask = [];
roiMasks=[];

for DepthToLoad=1:nDepthsTotal;
    fprintf(['Extracting Depth: ' num2str(DepthToLoad) '.\n'])
    load(fullfile(path,[loadList{DepthToLoad} '.mat']),'dat');
    cellID = [dat.stat(:).iscell];
    
    yoffset = dat.ops.yrange(1);
    xoffset = dat.ops.xrange(1);
    
    %Motion correct
    prev = s2pEpoch-1;
    Tstart = 1;
    while prev>0;
        Tstart = Tstart + size(dat.Fcell{prev},2);
        prev=prev-1;
    end
    
    
    T{DepthToLoad}= dat.ops.DS(Tstart:end,:);
    
    temp = dat.Fcell{s2pEpoch};
    if ~skipSP
        tempSP =dat.sp{s2pEpoch};
    end
    
    numCell = sum(cellID'==1);
    
    
    Fdata = temp(find(cellID'==1),:);
    if ~skipSP
        SPdata = tempSP(find(cellID'==1),:);
    end;
    Fnp = dat.FcellNeu{s2pEpoch}(find(cellID'==1),:);
    npc = [dat.stat(find(cellID'==1)).neuropilCoefficient];
    FdataNoNP = Fdata;
    Fdata= Fdata-Fnp.*npc';
    
    tempCOM = extractfield(dat.stat(cellID'==1),'med'); %get the center of each ROI
    tempCOM = reshape(tempCOM,2,numCell)'+[yoffset xoffset];
    
    %     Mask{DepthToLoad} = zeros([512 512]);
    
    ny=numel(dat.ops.yrange);
    nx=numel(dat.ops.xrange);
    
    thisPlane=zeros(ny,nx);
    thisPlane=thisPlane(:);
    cellvce=find(cellID'==1);
    for i = 1:numel(cellvce)
        %         if cellID(i)
        %             Mask{DepthToLoad}(dat.stat(i).ypix+yoffset,dat.stat(i).xpix+xoffset)=1;
        thisPlane(dat.stat(cellvce(i)).ipix,i)=1;
        %         end
    end
    
    theROIs=reshape(thisPlane,ny,nx,i);
    theROIImg=zeros(512,512,i);
    theROIImg(dat.ops.yrange, dat.ops.xrange,:)=theROIs;
    roiMasks = cat(3,roiMasks, theROIImg);
    
    %compile fluorescence data
    sz=size(Fdata);
    
    FdataNew=[];FdataNoNPNew=[];
    for i=1:numel(frameCount)
        FdataNew{i} = Fdata(:,sum(frameCount(1:i-1))+1 : min(sum(frameCount(1:i)),sz(2)) );
        FdataNoNPNew{i} = FdataNoNP(:,sum(frameCount(1:i-1))+1 : min(sum(frameCount(1:i)),sz(2)) );
        FNPNew{i} = Fnp(:,sum(frameCount(1:i-1))+1 : min(sum(frameCount(1:i)),sz(2)) );

    end
    
    temp=[];
    temp=cellfun(@(x) x(:,1:numFrames),FdataNew,'UniformOutput',false);
    Fdata = cell2mat(temp);
    Fdata = reshape(Fdata,sz(1),numFrames,numel(frameCount));
    
    temp=[];
    temp=cellfun(@(x) x(:,1:numFrames),FdataNoNPNew,'UniformOutput',false);
    FdataNoNP = cell2mat(temp);
    FdataNoNP = reshape(FdataNoNP,sz(1),numFrames,numel(frameCount));
    
    temp=[];
    temp=cellfun(@(x) x(:,1:numFrames),FNPNew,'UniformOutput',false);
    Fnp = cell2mat(temp);
    Fnp = reshape(Fnp,sz(1),numFrames,numel(frameCount));
    
    allDataNoNP = cat(1,allDataNoNP,FdataNoNP);
    allData = cat(1,allData,Fdata);
    allNP = cat(1,allNP,Fnp);
    allNPC = [allNPC npc];
    
    %compile Spiking Data
    if ~skipSP
        SPdataNew=[];
        for i=1:numel(frameCount)
            SPdataNew{i} = SPdata(:,sum(frameCount(1:i-1))+1 : min(sum(frameCount(1:i)),sz(2)) );
        end
        
        temp=[];
        temp=cellfun(@(x) x(:,1:numFrames),SPdataNew,'UniformOutput',false);
        SPdata = cell2mat(temp);
        SPdata = reshape(SPdata,sz(1),numFrames,numel(frameCount));
        
        allSPData = cat(1,allSPData,SPdata);
    end
    
    allCoM = cat(1,allCoM,tempCOM);
    allDepth = cat(1,allDepth,repmat(DepthToLoad,[numCell 1]));
end

backupData=allData;
if ~skipSP
    backupSPData=allSPData;
else
    backupSPData=[];
end
%Dimension order: ROI, Frame, Trial
sz = size(allData);


if UseNPCorrection
    disp('Processing With Neuropil Correction')
    dataToUse = allData;
else
        disp('Processing Without Neuropil Correction')
    dataToUse = allDataNoNP;
end


B = round(BaseLinePeriod/1000*FR);
BL = mean(mean(dataToUse(:,1:B,:),3),2);

% tempData = reshape(allData,[size(allData,1) size(allData,2)*size(allData,3)] );
% BL = prctile(tempData,3,2);

zProcess = 1;

allDataCorrect = dataToUse-min(BL);
% if zProcess == 0
%     BL= BL-min(BL);
%     DFFdata = (allDataCorrect-BL)./BL;
% elseif zProcess == 1 %zscore
%     DFFdata =zscore(allDataCorrect,0,2);
% elseif zPreocess ==2 %zscore dFF
%     BL= BL-min(BL);
% DFFdata = (allDataCorrect-BL)./BL;
% DFFdata =zscore(DFFdata,0,2);
% %
% end

BL= BL-min(BL);
DFFdata = (allDataCorrect-BL)./BL;

unroll = reshape(allDataCorrect,[sz(1) sz(2)*sz(3)]);
zData = zscore(unroll,0,2); %full recording zscore
zData = reshape(zData,sz);

zData2 =zscore(allDataCorrect,0,2); %trialwise zscore

if ~skipSP
    unroll = reshape(allSPData,[sz(1) sz(2)*sz(3)]);
    zSPData = zscore(unroll,0,2);
    zSPData = reshape(zSPData,sz);
    
    zSPData2 = zscore(allSPData,0,2);
    
     spBL = [];%mean(mean(allSPData(:,1:round(FR*2),:),3),2);
    
    if ~all(spBL>=0)
        disp('ERROR!! Negative spike counts')
    end
    
    miSPData = (allSPData-BL)./(allSPData+BL);
else
    zSPData=[];
    zSPData2=[];
    spBL=[];
    miSPData=[];
end


disp('Extracting Trial Information')

[stimID, uniqueStims, uStimCount] = stimReader(ExpStruct,DAQepoch);
% [runVector] = runReader(ExpStruct,DAQepoch,FR);
try
    [runVector] = runReader(ExpStruct,DAQepoch,FR);
catch
    
    ExpStruct.digitalSweeps=DIGITSWEEP;
    
    %     [runVector] = runReader2(ExpStruct,DAQepoch,FR,5); %normal
    [runVector] = runReader2(ExpStruct,DAQepoch,FR,3); %Satsuma Rig
    
end

numCells = sz(1);
%     imagesc(squeeze(Fdata(1,:,stimID==14))')
%%Determine Artifacts

disp('Determine Artifacts')
ArtifactSizeLeft = 130; %115 traditional
ArtifactSizeRight = 100;

ROIinArtifact = allCoM(:,2)<ArtifactSizeLeft | allCoM(:,2)>511-ArtifactSizeRight;
ROIs = 1:numel(ROIinArtifact);

excludeCells=ones(size(ROIs));

%%Exclude cells if nescessary
CellstToExclude =[];[9 8 6 ];
excludeCells = ~ismember(ROIs,CellstToExclude);
excludeCells = excludeCells' & ROIinArtifact;


%%Limit T to relevant frames for use with continuous imaging
disp('Reprocessing Transform Matrix')
framesUsed = 1:numFrames;
for i= 2:numel(frameCount)
    strtFrame = sum(frameCount(1:i-1));
    framesUsed = [framesUsed strtFrame:strtFrame+numFrames-1];
end

for k = 1:numel(T)
    Tused{k}= T{k}(round(framesUsed),:);
end