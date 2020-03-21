%% Load one Out

%% Identify Red Cells from another image
[redImageFN redImagePath] = uigetfile;
% rawRedImg = bigread3([redImagePath redImageFN]);
rawRedImg = ScanImageTiffReader([redImagePath redImageFN]).data;
rawRedImg=rawRedImg(:,:,2:2:end);

% exampleGreenImage = ScanImageTiffReader(localFiles{1}).data;
% exampleGreenImage = exampleGreenImage(:,:,1:2:end);
%% Manually click Red Cells
nDepthsTotal = max(out.exp.allDepth);

    RedCellLocations=[];

for i = 1:nDepthsTotal
%     subplot(1,3,i)
figure(7);clf
imToUse =mean(rawRedImg(:,:,i:nDepthsTotal:end),3);
    imagesc(imToUse);
    title(['Plane ' num2str(i) '. Click on every ''Red'' Cell. click outside the block to moveon']);
    
    caxis([prctile(imToUse(:),0.5) prctile(imToUse(:), 99.5)]);
    hold on

    thisLocation =[0 0];
    while all(thisLocation>=0) && all(thisLocation<=512)
        thisLocation = ginput(1);
        plot(thisLocation(1),thisLocation(2),'o','color',rgb('red'))
        RedCellLocations = cat(1,RedCellLocations,[thisLocation i]);
    end
        
    
%     plot(allCoM(allDepth==i,1)-offsets(1),allCoM(allDepth==i,2)-offsets(2),'o','color',rgb('red'))
%     plot(stimCoM(stimDepth==i,1),stimCoM(stimDepth==i,2),'.')
    xlim([0 512])
    ylim([0 512])
end
clf
%% retrieve saved varialbes
offsets = out.info.offsets;
allDepth = out.exp.allDepth;
allCoM = out.exp.allCoM;
cellList = 1:numel(allDepth); 
numCells = numel(allDepth); 
HoloTargets = out.exp.holoTargets;
additionalOffsets = zeros([nDepthsTotal 2]);

%% optional additional offsets.
%offset is just a single average number, but sometimes you can tell that
%one plane is more offset then the rest ad that now

additionalOffsets(1,:)=[1 2];
additionalOffsets(2,:)=[0 1];
additionalOffsets(3,:)=[1 1];

%% match red cell locations
figure(7);clf
redDepth = RedCellLocations(:,3); 
RedCellLoc = nan(size(RedCellLocations));
for i = 1:nDepthsTotal
    offsetToUse =offsets + additionalOffsets(i,:); 
        subplot(1,3,i)
        plot(allCoM(allDepth==i,1)-offsetToUse(1),allCoM(allDepth==i,2)-offsetToUse(2),'o','color',rgb('blue'))
hold on
plot(RedCellLocations(redDepth==i,1),RedCellLocations(redDepth==i,2),'o','color',rgb('red'))
    xlim([0 512])
    ylim([0 512])
    
    RedCellLoc(redDepth==i,:) = RedCellLocations(redDepth==i,:)-[offsetToUse 0];

end

numRed=size(RedCellLoc,1);

RedMapping =zeros([numRed 3]);
for s = 1:numRed
    
    d = RedCellLoc(s,3);
    RedMapping(s,1) = d;
    
    a = allCoM;
    a(allDepth~=d,:) = nan;
    
    b = RedCellLoc(s,1:2);
    [RedMapping(s,3), RedMapping(s,2)] = min ( sqrt(sum((a-b).^2,2)) );
end


targetDistanceThresholdRed =15; %has been 15 
RedCells = RedMapping(:,2); 
RedCells(RedMapping(:,3)>targetDistanceThresholdRed)=nan;
isRed = ismember(cellList,RedCells);

% disp('Calculating RedValues')
% redVal = zeros([numCells 1]);
% for i=1:numel(RedCells)
%     if ~isnan(RedCells(i))
%         r=RedCells(i);
%         d = allDepth(r);
%         im = mean(rawRedImg(:,:,d:nDepthsTotal:end),3);
%         redVal(r) = mean(unroll(im(boolean(roiMasks(:,:,r)))));
%     end
% end

disp(['Found ' num2str(sum(isRed)) ' Red Cells. Out of ' num2str(numel(RedCells)) ' Targets clicked.'])
targetedRedCells = find(ismember(cellList,[HoloTargets{:}]) & isRed);
problemHolos = cellfun(@(x) any(ismember(x,targetedRedCells)),HoloTargets,'uniformoutput',1);

disp([num2str(numel(targetedRedCells)) ' Red Cells were used in Holos.'])

%% save to out
red.RedCellLoc = RedCellLoc;
red.isRed = isRed;
% red.redVal = redVal; 
red.RedCells = RedCells;

out.red=red;

info = out.info;

save(['U:\ioldenburg\outputdata1\' info.date '_' info.mouse '_outfile'], 'out')
disp('saved')

