function   [localFiles, frameCount, MD, FR, minNumFrames, numFrames, numCells,...
    numFiles, allData, allSPData, allCoM, allDepth, allNP, allNPC,...
    roiMasks, T, BL, allDataCorrect, zData, zSPData, backupData,...
    DFFdata, backupSPData, ExpStruct, CellstToExclude, allDataNoNP,...
    stimID, uniqueStims, uStimCount, runVector, ArtifactSizeLeft,...
    ArtifactSizeRight, ROIinArtifact, excludeCells, framesUsed, Tused, zData2]...
    = loadCaDataPY(path,DAQepoch,s2pEpoch,nDepthsTotal,ExpStruct,DIGITSWEEP,BaseLinePeriod,numColors)


UseNPCorrection=1;

dat = load(path{1});
epoch_starts = find(dat.ops.first_tiffs);

% this will get the number of frames per tiff
try
    frameCount = dat.ops.frames_per_file(epoch_starts(s2pEpoch):epoch_starts(s2pEpoch+1));
catch
    frameCount = dat.ops.frames_per_file(epoch_starts(s2pEpoch):end);
end



frameCount(frameCount==1)=[];

% use to get metadata and framerate
localFile = dat.ops.filelist(1,:);

if ~exist('im')
    im = imfinfo(localFile);
end

MD = parseSI5Header(im(1).Software);
FR = MD.SI.hRoiManager.scanVolumeRate;
minNumFrames=min(frameCount);
numFrames =floor(minNumFrames);
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

for DepthToLoad=1:nDepthsTotal
    fprintf(['Extracting Depth: ' num2str(DepthToLoad) '.\n'])
    dat = load(fullfile(path{DepthToLoad}));
    
    % split F, Fneu, and spikes into epochs
    dat.Fcell = mat2cell(dat.F, size(dat.F,1), dat.ops.frames_per_folder);
    dat.FcellNeu = mat2cell(dat.Fneu, size(dat.Fneu,1), dat.ops.frames_per_folder);
    dat.sp =  mat2cell(dat.spks, size(dat.spks,1), dat.ops.frames_per_folder);
    
    cellID = dat.iscell(:,1)';
    
    dat.ops.yrange = [dat.ops.yrange(1):dat.ops.yrange(2)];
    dat.ops.xrange = [dat.ops.xrange(1):dat.ops.xrange(2)];
    
    yoffset = dat.ops.yrange(1);
    xoffset = dat.ops.xrange(1);
    
    %Motion correct
    prev = s2pEpoch-1;
    Tstart = 1;
    while prev>0
        Tstart = Tstart + size(dat.Fcell{prev},2);
        prev=prev-1;
    end
    
    % this might be backwards?? x-y or y-z?
    dat.ops.DS = [dat.ops.xoff; dat.ops.yoff]';
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
%     npc = [dat.stat(find(cellID'==1)).neuropilCoefficient];
    npc = dat.ops.neucoeff;
    FdataNoNP = Fdata;
    Fdata= Fdata-Fnp.*npc';
    
    dat.stat = cell2mat(dat.stat);
    
    tempCOM = extractfield(dat.stat(cellID'==1),'med'); %get the center of each ROI
    tempCOM = reshape(tempCOM,2,numCell)'+double([yoffset xoffset]);
    
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
        disp(['cell = ' num2str(i) ', size= ' num2str(size(thisPlane))])
        %         end
    end
    
    %% reshape throwing an error...
    % numel(thisPlane) != ny*nx*i
    
    theROIs=reshape(thisPlane,ny,nx,i);
    theROIImg=zeros(512,512,i);
    theROIImg(dat.ops.yrange,dat.ops.xrange,:)=theROIs;
    roiMasks = cat(3,roiMasks, theROIImg);
  %%  
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