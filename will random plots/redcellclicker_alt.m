% [loadList, loadPath ]= uigetfile('E:\100spikes-results\outfiles-master','MultiSelect','on');]

%%

numExps = numel(loadList);
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

out = All.out;
out.info
nDepthsTotal=3;

%% Identify Red Cells from another image (INTERNEURON)

%get the 830nm Image
[redImageFN redImagePath] = uigetfile('*.tif');
rawRedImg = ScanImageTiffReader([redImagePath redImageFN]).data;
rawRedImg=rawRedImg(:,:,2:2:end);

%% Get 1020 image too

[red1020FN redImagePath] = uigetfile([redImagePath '*.tif']);
rawRed1020Img = ScanImageTiffReader([redImagePath red1020FN]).data;
rawRed1020Img=rawRed1020Img(:,:,2:2:end);

%%
%% Motion Correct and check each

clear IMTOUSE
for i = 1:nDepthsTotal
    
    imToUse = rawRedImg(:,:,i:nDepthsTotal:end);
    imToUse = mean(imToUse,3);
    imToUse = imToUse - min(imToUse, [], 'all');
    
    imToUse2 = rawRed1020Img(:,:,i:nDepthsTotal:end);
    imToUse2 = mean(imToUse2,3);
    imToUse2 = imToUse2 - min(imToUse2, [], 'all');
    
    IMTOUSE1{i} =  int16(imToUse);
    IMTOUSE2{i} =  int16(imToUse2);
    
end

%% Manually Click Red Cells Using merge
nDepthsTotal = max(out.exp.allDepth);

RedCellLocations=[];

for i = 1:nDepthsTotal
    
    figure(7)
    clf
    
    imToUse = IMTOUSE1{i};
    imToUse2 = IMTOUSE2{i};
    
    IMG = zeros([512 512 3]);
    IMG(:,:,1) = imToUse2;
    IMG(:,:,2) = imToUse;
    IMG(:,:,3) = imToUse;
    image(IMG);
    
    title(['Plane ' num2str(i) '. Click on every blue Cell. click outside the block to moveon']);
    hold on
    
    thisLocation =[0 0];
    while all(thisLocation>=0) && all(thisLocation<=512)
        thisLocation = ginput(1);
        plot(thisLocation(1),thisLocation(2),'o','color',rgb('red'))
        RedCellLocations = cat(1,RedCellLocations,[thisLocation i]);
    end
    
    xlim([0 512])
    ylim([0 512])
end
clf
