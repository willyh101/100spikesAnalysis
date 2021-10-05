%% Start Here:
% Hi Lamiae, this is my code for making movies. It comes at the end of one
% of my analysis files so i'm not sure if all the file formats will be
% helpful. It very much believes you to be imaging 3 and only 3 planes. I
% don't think i've ever run it with 1, so i'm not sure how it will work.

%Its best to tryt to read through the whole code before running. I've tried
%to anotate but i'm not sure whats ian specific vs universal.

%Variables You will need:
%lowMotionTrials- binary list of trials in which there was less than some
%threshold of movement (used for data exclusion)
%highRunTrials - binary list of trials in which the mouse was running above
%threshold (used for data exclusion)
%stimID - a trialwise list of each stimulation type (so 1 would be no stim, 2
%the first holo, etc)
%visID - a trialwise list of each visual stimulus given each trial (1 is no
%vis, 2 is first grating or contrast, etc)
%tTest - is a cell array for each depth of the per frame movement (its median
%subtracted, and only includes the frames that are used in the rest of the
%analysis, aka numFrames)
%numFrames - i crop all my movies to a minimal number of frames
%strt - start time of holostim in ms
%dura - duration of holostim in ms
%freq - if more than one holostim per trial the inter holo interval in hz
%stimParam.roi - a cell array with each cell giving an ID of the hologram
%used in that stimID. so trial 1 displays the 5th hologram that
%contains 3 targets; therefore stimID(1)=5 so stimParam.roi{5} =
%ROIID; roiTargets{ROIID} = [targ1 targ2 targ3]; holoTargets{ROIID} =
%[matchedCell1 matchedCell2 matchedCell3];
% roisTargets - cell array of each holograms targets shot (numbers
% correspond to ROIs which are located as in stimCoM;
% holoTargets - as with roisTargets but translated to matched cells, ie
% indexed by allCoM
%StimCoM and stimDepths - the location and depth of each stimulus (in SI
%pixels
%AllCoM and allDepths - as above for every detected cell.
%offsetts - the median movement translation to get everything to line up
%right
%ArtifactSizeLeft and ArtifactSizeRight - size of your artifacts usually
%100 and 100


% I think thats it. Good Luck...

%%
oldLocalFiles = localFiles;
%%
for i=1:numel(localFiles);
%     localFiles{i} = ['Z:\holography\Data\ian\180523\' oldLocalFiles{i}(57:end)];
%     localFiles{i} = ['Z:\holography\Data\ian\180522\' oldLocalFiles{i}(59:end)];
    localFiles{i} = ['Z:\holography\Data\ian\180522\' oldLocalFiles{i}(59:end)];
localFiles{i} = ['T:\Will\20201103\w29_1\'  oldLocalFiles{i}(42:end)];
end
%% Find and Motion Correct Raw data (very Slow)

d1Mov=[];
d2Mov=[];
d3Mov=[];


trialsToUse = lowMotionTrials & ~highRunTrials;

uniqueStims = unique(stimID);
uniqueVisStim = unique(visID);
% trialsToUse = all(trialsToUse);

i=0;
for idx = 1:numel(uniqueStims)
    for k = 1:numel(uniqueVisStim)
        i=i+1;
        disp(['Cycle ' num2str(i) ' of ' num2str(numel(uniqueStims)*numel(uniqueVisStim)) ]);
        ID = uniqueStims(idx);
        vID = uniqueVisStim(k);
        criteria = trialsToUse & stimID==ID & visID ==vID;
        log(i,:)=[ID vID];
        % criteria = ones([1 numFiles]);
        
        if sum(criteria)>0
            [allMovs]= meanMovieMakerN(...
                localFiles,...
                criteria,...
                'g',...
                cellfun(@(x) -round(x),tTest,'uniformoutput',0),...%cellfun(@(x) -round(x),Tused,'uniformoutput',0),...
                numFrames,...
                nDepthsTotal);
            for j = 1:nDepthsTotal
                eval(['d' num2str(j) 'Mov{i}= allMovs{' num2str(j) '};'])
            end
            
        else
            for j = 1:nDepthsTotal
                eval(['d' num2str(j) 'Mov{i}= [];'])
            end
            %         d1Mov{i}= [];
            %         d2Mov{i}=[];
            %         d3Mov{i}=[];
        end
    end
end

%% Compute PreZScored Movies (a bit slow)

mov1 = cat(3,d1Mov{:});
mov2 = cat(3,d2Mov{:});
mov3 = cat(3,d3Mov{:});
mov4 = cat(3,d4Mov{:});


disp('Computing ZScore 1')
zMov1 = zscore(mov1,0,3);
disp('Computing ZScore 2')
zMov2 = zscore(mov2,0,3);
disp('Computing ZScore 3')
zMov3 = zscore(mov3,0,3);
disp('Computing ZScore 4')
zMov4 = zscore(mov4,0,3);
disp('done')

%% Play Mov

saveVid = 1; %save or do not save, there is no try
zScoreIt=0; %trialwise zscore, faster but less accurate
preZScore =1; %Can Take a long time, not to be used with zScoreIt
gausfilt = 1; %its helpful to gausian fliter it to reduce noise
gausfiltVal = 0.75; % How much to gaussian filter typically 0.75 sometimes 1.25
baselineSubtract = 1; %subtract a prestimulus period
cellsOnly = 0; % only the outlines of cells
blankEdges = 1; % put white bars over where the artifacts will be
trimEdges = 1; %simply cut off where the artifacts will be
rotateIM = 1; %rotate the image 90degrees so that it looks landscape rather than portrait after edge trimming
useStimCoM=1; %use the location of stimulated holograms, rather than their matched cells to place circles
superimposeIMs = 0; %collapse all planes into a single one
maxFrameDisplay = inf; 18; %Cut off Movie early if you want, default inf
oneFrameAtATime = 0; % this will pause after each frame if you're trying to pull out a single image

stimToUseList=[2 5 4];%[2 5 4];[ 8 26  ];% %[1 4 10 2 5 11];%[1 4 10 2 5 11];% 1:15% 2  4 5 8];%[7 13] ;[1 2 3]; %[4 7 10 13] ;[7 8 9]; %[4 5 6]; %[1 2 3];% The StimIDs that you want to cycle through

aF = 1; %average Factor; How many frames do you want average. default 1

useVisStimMarker = 1; %display a marker when the visual stimulus is on.
useVisStimMarkerText =1;% display text describing that visual stimulus
% stimMarkerText = {'Movie Broken' '0 Contrast' '0.01 Contrast' '0.04 Contrast' '0.1 Contrast' '0.4 Contrast' '1 Contrast'};
stimMarkerText = { 'No Grating' '90\circ Grating' '135\circ Grating'};
stimMarkerText = {'No Vis' '0\circ' '45\circ' '90\circ' '135\circ' '180\circ' '225\circ' '270\circ' '315\circ'};

useHoloStimMarkerText = 1; %write the holo number onto it

visMarkerSize = 25;
visMarkerColor =  rgb('Magenta');


stimstyle = 2;  % if you're not doing any visual stimulus stimstyle 1 and 2 should be the same, otherwise you probably want to use 2
numHolosPerTrial=1; %if a given sweep (trial type) has more than one hologram set that here. This was last used with only 1 type of hologram and might not work anymore with multiple.

extendStimulus = 5; %frames to extend stim by
extendColor = rgb('violet');%rgb('grey');

offPlaneStim = 0; %draw a circle in the off plane locations
offPlaneStimColor = rgb('darkgray');

drawscalebar = 1; %draw a scalebar
scaleBarColor = rgb('LightSlateGray');% rgb('White')

playSpeed = 2;%2.5; %Playback faster? only matters for watching not for saving. multiple of FrameRate

displayTitle = 1; %display title of stimulus (stimID)

depthString = {'-390\mum' '-350\mum' '-193\mum' '-150\mum'};%-60\mum' '-30\mum' '0\mum'}; % Title of each plane
% depthString = {'-188\mum' '-158\mum' '-128\mum'};%-60\mum' '-30\mum' '0\mum'}; % Title of each plane
clim =[0 2];[-4 4]; % default colormap
asymetricCLIM = 0; % sometimes you want a colormap that isn't smooth, but has different movement <0 and >0 this allows that

labelLayers = 1; %cludgey to label the first two plots L4 and the second L23
%set subplots
co = 2;
ro = 2;2
nDepthsTotal = 4;4;
% clim = [-100 200];[0 2.5];


circleTunedCells = 0; %circle tuned cells
circleCriteria = find(pVisR<0.05 & prefVR==8);%find(pVisR<0.05 & prefVR==2);
circleDiameter = 8;
circleColor = rgb('cyan');
circleLineWidth =1;

% clim = [0 4];%4];%2.5];
%  clim = [-3 3];
if superimposeIMs
    clim=clim;
end

if asymetricCLIM
    ratio= abs(clim(2)/clim(1));
    temp = rdbu; %Set the Colormap here
    szCM = size(temp);
    temp1=temp(1:ceil(szCM(1)/2),:);
    temp2= temp(ceil(szCM(1)/2)+1:szCM(1),:);
    nCM = size(temp2,1);
    temp3 = [interp1(temp2(:,1),1:1/ratio:nCM); interp1(temp2(:,2),1:1/ratio:nCM); interp1(temp2(:,3),1:1/ratio:nCM)]';
    temp4 = [temp1; temp3];
    % maptouse = b2r(-3.5,3.5,[temp(1,:); [1 1 1]; temp(end,:)]);
    maptouse = temp4;
else
    maptouse = viridis; rdbu; flipud(AdvancedColormap('Royal')); viridis; rdbu; viridis; %Set the colormap here rdbu is what i use most often
end

%how to make the stimulus marker look
stimColor = rgb('Magenta');rgb('Magenta');%rgb('MistyRose') ;
stimLineWidth = 1.5;
stimDiameter = 10;

%save video parameters
if saveVid
    savePath = 'C:\Users\ian\Documents\DATA\AnalysisOutput';
    %         savePath = 'C:\Users\SabatiniLab\Dropbox\Adesnik\Data To Play With';
    
    name = '201124_w29_1_L4toL23';
    
    vid = VideoWriter(fullfile(savePath,[name '.avi']));
    vid.FrameRate =FR*playSpeed;  % Default 30
    vid.Quality = 100;    % Default 75
    
    open(vid)
end

movFrames = size(d1Mov{stimToUseList(1)},3);


clear h1
for M=1:numel(stimToUseList)
    stimToUse = stimToUseList(M);
    
    %all the beautification code
    if preZScore
        ind = stimToUse; %this might need to change depending on how zMov was built
        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) ' = zMov' num2str(j) '(:,:,(ind-1)*movFrames+1:ind*movFrames);']);
        end
    else
        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) ' = d' num2str(j) 'Mov{stimToUse};']);
        end
    end
    
    if zScoreIt

        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) '=  zscore(mov' num2str(j) ',3,0,3);']);
        end
    else
        %         clim = [-25 25];
    end
    
    if gausfilt

        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) '=  imgaussfilt(mov' num2str(j) ',gausfiltVal);']);
        end
    end
    
    strtFrame = round(strt/1000*FR)-1;
    if baselineSubtract

        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) '= mov' num2str(j) ' - imgaussfilt(mean(mov' num2str(j) '(:,:,1:strtFrame),3),2.5);']);
        end
    end
    
    if cellsOnly
        crit = ones([numel(ROIinArtifact) 1]);%~ROIinArtifact;

        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) '= mov' num2str(j) '.* any(roiMasks(:,:,allDepth(1:numCells)==' num2str(j) ' & crit),3);']);
        end
    end
    
    blankMov = zeros(size(mov1));
    if blankEdges

        for j = 1:nDepthsTotal
            temp = blankMov;
            eval(['temp(:,ArtifactSizeLeft:511-ArtifactSizeRight,:)= mov' num2str(j) '(:,ArtifactSizeLeft:511-ArtifactSizeRight,:);']);
            eval(['mov' num2str(j) ' = temp;']);
        end
        
    end
    
    
    %section to determine when and where to display a stim
    when =[];% zeros([1 size(sMeanMov,3)]);
    whenLong = zeros([1 size(mov1,3)]);
    whenLongVis = zeros(size(whenLong));
    
    for i=1:numHolosPerTrial
        
        stimTime= strt/1000 + (i-1)/freq; %onset of stim in s
        stimLen = dura/1000; %duration of stim in s
        
        strtStim = round((stimTime*FR));
        endStim = round((stimTime*FR)+stimLen*FR);
        
        when = [when strtStim:endStim];
        % whenLong(strtStim:endStim)=i;
        if stimstyle==1
            whenLong(strtStim:endStim)=i;%stimParam.roi(stimToUse); %for one stimulus desing
        elseif stimstyle ==2
            stimThisMov = find(log(stimToUse,1)==uniqueStims);
            whenLong(strtStim:endStim)=    stimParam.roi{stimThisMov};
            visThisMov = find(log(stimToUse,2)==uniqueVisStim);
            whenLongVis(strtStim:endStim) = visThisMov;
        end
    end
    
    whenLong2 = zeros(size(whenLong));
    if extendStimulus>0
        ns = find(whenLong~=0);
        
        for i =1:numel(ns);
            n = ns(i);
            stim = whenLong(n);
            whenLong2(n+1:n+extendStimulus)=stim;
        end
    end
    
    if useStimCoM
        if size(allCoM,1)==numCells
            allCoMBackup = allCoM;
            allCoM2 = cat(1,allCoM,stimCoM);
            allDepthBackup = allDepth;
            allDepth2 = cat(1,allDepth,stimDepth);
        elseif size(allCoM,1)~=numCells+size(stimCoM,1);
            errordlg('Error in the size of all CoM');
        else
            allCoM2 = allCoM;
            allDepth2 = allDepth;
        end
        
        tempTargets = cellfun(@(x) x+numCells, roisTargets,'UniformOutput',0); %HolosTargeted is the variable containing the stim targets
    else
        
        allCoM2 = allCoM;
        allDepth2 = allDepth;
        
        tempTargets = HoloTargets;
        for i=1:numel(HoloTargets)
            tempTargets{i}(isnan(tempTargets{i}))=[];
        end
    end
    allCoM2 = allCoM2-offsets;
    
    
    if rotateIM

        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) '= permute(mov' num2str(j) ',[2 1 3]);']);
        end
        allCoM2 = fliplr(allCoM2);
    end
    
    if trimEdges

        for j = 1:nDepthsTotal
            eval(['mov' num2str(j) '= mov' num2str(j) '(ArtifactSizeLeft:511-ArtifactSizeRight,:,:);']);
        end
        allCoM2 = allCoM2-[ArtifactSizeLeft 0];
        allCoM2(allCoM2(:,1)<0,:)=NaN;
    end
    
    if superimposeIMs
        temp = mov1;
        for j=2:nDepthsTotal
            eval(['temp = cat(4,temp,mov' num2str(j) ');']);
        end
        mov1 = mean(temp,4)*nDepthsTotal;
        
    end
    
    
    F = figure(115);clf
    imagesc(mov1(:,:,1));
    set(gcf,'color','w')
    hold on
    axis tight
    axis square
    axis off
    h = colorbar;
    
    colormap(maptouse)
    
    if superimposeIMs
        %     F.Position = [-1154 1006 986 504];
    else
        
    end
    
    %
    for  i =  aF+1:min(size(mov1,3),maxFrameDisplay);
        clf
        if superimposeIMs
            imagesc(mean(mov1(:,:,i-aF:i),3));
            caxis(clim);
            axis equal
            axis off
            h = colorbar;
            colorbarPos = h.Position;
            ylabel(h,'\Delta Z-Score Fluorescence');
            set(h.Label,'rotation',-90);
            
            if trimEdges
                h.Position=colorbarPos+[+.05 -.02 0 .04];
                h.Label.Position = [3. mean(clim) 0];
            else
                h.Position=colorbarPos+[-.14 -.02 0 .05];
                h.Label.Position = [3 mean(clim) 0];
            end
            
            h.Label.FontSize = 14;
            
            if whenLong(i)~=0
                cells = tempTargets{whenLong(i)};
                for k=1:numel(cells)
                    cell = cells(k);
                    
                    temp = circle(allCoM2(cell,2),allCoM2(cell,1),stimDiameter,stimColor);
                    temp.LineWidth = stimLineWidth;
                    temp.LineStyle = ':';
                    hold on
                end
            end
            
            if useVisStimMarker & whenLongVis(i)~=0
                visMarker = line(20,20);%,stimDiameter,stimColor);
                visMarker.LineStyle = 'none';
                visMarker.Color = visMarkerColor;
                visMarker.MarkerFaceColor = visMarkerColor;
                visMarker.MarkerSize=visMarkerSize;
                visMarker.Marker = 's';
            end
            
            
            if extendStimulus>0
                if whenLong2(i)~=0
                    cells = tempTargets{whenLong2(i)};
                    for k=1:numel(cells)
                        cell = cells(k);
                        
                        temp = circle(allCoM2(cell,2),allCoM2(cell,1),stimDiameter,extendColor);
                        temp.LineWidth = stimLineWidth;
                        temp.LineStyle = ':';
                        hold on
                    end
                end
            end
            
            if drawscalebar
                if trimEdges
                    r = rectangle('Position',[10 270 128 2]);
                    r.LineWidth =2;
                    r.EdgeColor = scaleBarColor;
                    r.FaceColor = scaleBarColor;
                    
                    t = text(128/2,270-10,'200\mum');
                    t.Color=scaleBarColor;
                    t.FontSize=14;
                else
                    r = rectangle('Position',[10 500 128 2]);
                    r.LineWidth =2;
                    r.EdgeColor = scaleBarColor;
                    r.FaceColor = scaleBarColor;
                    
                    t = text(128/2-10,500-15,'200\mum');
                    t.Color=scaleBarColor;
                    t.FontSize=14;
                end
            end
            if displayTitle
                if stimToUse==1
                    title('No Stimulation')
                else
                    title(['Stimulus : ' num2str(stimToUse)])
                end
            end
        else %don't superimpose
            
            for j=1:nDepthsTotal
                
                h1{j} = subplot(ro,co,j);
                eval(['imagesc(mean(mov' num2str(j) '(:,:,i-aF:i),3));']);
                caxis(clim);
                if ~trimEdges
                    axis square
                    title(depthString{j},'Units','normalized','Position',[0.5 1 1])
                else
                    title(depthString{j},'Units','normalized','Position',[0.5 1 1])
                    axis equal
                end
                axis off
                HoloTxt=[];
                if whenLong(i)~=0
                    cells = tempTargets{whenLong(i)};
                    for k=1:numel(cells)
                        cell = cells(k);
                        
                        if allDepth2(cell)==j
                            temp = circle(allCoM2(cell,2),allCoM2(cell,1),stimDiameter,stimColor);
                            temp.LineWidth = stimLineWidth;
                            temp.LineStyle = ':';
                            hold on
                        elseif offPlaneStim & (allDepth2(cell)==j-1 | allDepth2(cell)==j+1)
                            temp = circle(allCoM2(cell,2),allCoM2(cell,1),stimDiameter,offPlaneStimColor);
                            temp.LineWidth = stimLineWidth;
                            temp.LineStyle = ':';
                            hold on
                            
                        end
                    end
                    
                    
                     if useHoloStimMarkerText & j==1 & whenLong(i)>0
                        HoloTxt= ['Hologram: ' num2str(whenLong(i))]; %HoloStimMarkerText{whenLong(i)};
                        if ~useVisStimMarkerText
                        x = get(gca);
                        oldText = x.Title.String;
                        x.Title.String={HoloTxt; oldText};
                        end
                    end
                end
                
                if useVisStimMarker & whenLongVis(i)~=0
                    visMarker = line(20,20);%,stimDiameter,stimColor);
                    visMarker.LineStyle = 'none';
                    visMarker.Color = visMarkerColor;
                    visMarker.MarkerFaceColor = visMarkerColor;
                    visMarker.MarkerSize=visMarkerSize;
                    visMarker.Marker = 's';
                    
                    if useVisStimMarkerText & j==1
                        txt= stimMarkerText{whenLongVis(i)};
                        if useHoloStimMarkerText & ~isempty(HoloTxt)
                            txt = [txt ' | ' HoloTxt];
                        end
                        x = get(gca);
                        oldText = x.Title.String;
                        x.Title.String={txt; oldText};
                    end
                end
                
                if extendStimulus>0
                    if whenLong2(i)~=0
                        cells = tempTargets{whenLong2(i)};
                        for k=1:numel(cells)
                            cell = cells(k);
                            
                            if allDepth2(cell)==j
                                temp = circle(allCoM2(cell,2),allCoM2(cell,1),stimDiameter,extendColor);
                                temp.LineWidth = stimLineWidth;
                                temp.LineStyle = ':';
                                hold on
                            elseif offPlaneStim & (allDepth2(cell)==j-1 | allDepth2(cell)==j+1)
                                temp = circle(allCoM2(cell,2),allCoM2(cell,1),stimDiameter,offPlaneStimColor);
                                temp.LineWidth = stimLineWidth;
                                temp.LineStyle = ':';
                                hold on
                            end
                        end
                    end
                end
                
                if circleTunedCells 
                     cells = circleCriteria;
                    for k=1:numel(cells)
                        cell = cells(k);
                        
                        if allDepth2(cell)==j
                            temp = circle(allCoM2(cell,2),allCoM2(cell,1),circleDiameter,circleColor);
                            temp.LineWidth = circleLineWidth;
                            temp.LineStyle = ':';
                            hold on 
                        end
                    end
                end
                
                if drawscalebar & j==1
                    if trimEdges
                        r = rectangle('Position',[10 270 128 2]);
                        r.LineWidth =2;
                        r.EdgeColor = scaleBarColor;
                        r.FaceColor = scaleBarColor;
                        
                        t = text(128/2-15,270-15,'200\mum');
                        t.Color=scaleBarColor;
                        t.FontSize=14;
                    else
                        r = rectangle('Position',[10 500 128 2]);
                        r.LineWidth =2;
                        r.EdgeColor = scaleBarColor;
                        r.FaceColor = scaleBarColor;
                        
                        t = text(128/2-15,500-15,'200\mum');
                        t.Color=scaleBarColor;
                        t.FontSize=14;
                    end
                end
                
                
            end
            origPos = get(h1{end},'Position');
            
            h = colorbar;
            colorbarPos = h.Position;
            ylabel(h,'\Delta Z-Score Fluorescence');
            set(h.Label,'rotation',-90);
            h.Label.FontSize = 11;
            
%                         h.Position= max(colorbarPos+[0.03 .125 0 -.25],0);

            h.Position= max(colorbarPos+[0.06 .00 0 0],0);
            h.Label.Position = [2.5 mean(clim) 0];
            set(h1{end},'position',origPos);
            axis equal
        end
        
        if labelLayers
            subplot(ro,co,1);
            ax = gca;
            ylabel(ax,'Layer 4');
            ax.YLabel.Visible = 'on';
            
                        subplot(ro,co,3);
            ax = gca;
            ylabel(ax,'Layer 2/3');
            ax.YLabel.Visible = 'on';
        end
        
        if saveVid
            frame = getframe(gcf);
            writeVideo(vid,frame);
        else
            pause(1/(FR*playSpeed))
        end
        
        if mod(i,25)==0
            pause(0.01);
            fprintf('Paused 1s\n')
        end
        if mod(i,100)==0
            figure(115);clf
            fprintf('reset image\n')
        end
        
        fprintf(['Frame: ' num2str(i) '\n']);
        
        if oneFrameAtATime
            input('PressEnter')
        end
    end
    
    
    
end

if saveVid
    close(vid)
end

% clear zMov1 zMov2 zMov3

disp('done')
