%% This section Designed to Work with single experiment files

%% Orientation Curves

data=zdfData; dfData; zData;
c=0;
win =  round([visStart*FR (visStop)*FR]) + [2 2];%[8 17];%
visR=[];visS=[];visC=[];
clear X
for i = 1:numel(uniqueVisStim)
    v = uniqueVisStim(i);
    if v~=0
        for k = 1;%:numel(uniqueStims)
            c=c+1;
            s=uniqueStims(k);
            x = data(:,:,visID==v & stimID == s);
            X{i,k}= x;
%             
%             b2 = mean(x(:,1:3,:),2);
%             m2 = x-b2;
%             m = mean(m2,3);
%                         X{i,k}= m2;

            m = mean(x,3);
%            

            visR(c,:) = nanmean(m(:,win(1):win(2)),2);
            
             b = mean(m(:,1:3),2);
            m = m-b;
            visRB(c,:) = nanmean(m(:,win(1):win(2)),2);

            
            visS(c,:) = nanstd(squeeze(nanmean(x(:,win(1):win(2),:),2) ),[],2);
            visC(c,:) = size(x,3);
        end
    end
end

toDelete = [];
for i=1:numel(X)
    if isempty(X{i})
       toDelete = [toDelete i];
    end
end
X(toDelete)=[];

[maxVR prefVR] = max(visR,[],1);


pVisR=[];
for i = 1:numCells
        dat1 = cellfun(@(x) squeeze(mean(x(i,win(1):win(2),:))),X,'uniformoutput',0);
        tempDat = cat(1,dat1{:});
        tempGroup = [];
        for k =1:numel(dat1)
            tempGroup = [tempGroup ones([1 numel(dat1{k})])*k];
        end
        pVisR(i) = anovan(tempDat,tempGroup','display','off');%,'off');
        
        tempDat = cat(1,dat1{2:end});
        tempGroup = [];
        for k =2:numel(dat1)
            tempGroup = [tempGroup ones([1 numel(dat1{k})])*k];
        end
        pVisT(i) = anovan(tempDat,tempGroup','display','off');
end

alpha = 0.01;
disp(['Detected ' num2str(sum(pVisR<alpha)) ' Vis Responsive cells (p<' num2str(alpha) '). Approx ' num2str(sum(pVisR<alpha)/numCells*100,2) '%'])
disp(['Detected ' num2str(sum(pVisT<alpha)) ' Tuned cells (p<' num2str(alpha) '). Approx ' num2str(sum(pVisT<alpha)/numCells*100,2) '%'])

visRcells = pVisR<alpha;
visTcells = pVisT<alpha;
%% Compute Distance Matrix

Distance = [];

muPerPx = 800/512;

allLoc = [allCoM*muPerPx (allDepth-1)*30];
for i =1:numCells
    for k=1:numCells
        Distance(i,k) = sqrt(sum((allLoc(i,:)-allLoc(k,:)).^2));
        radialDistance(i,k) = sqrt(sum((allCoM(i,:)-allCoM(k,:)).^2)).*muPerPx;
        axialDistance(i,k) = sqrt(sum((allDepth(i)-allDepth(k)).^2));
    end
end
stimLoc = [stimCoM*muPerPx (stimDepth-1)*30];
%% Create Blast by Blast potentially ineligible Cells
thisPlaneTolerance = 12;10; %in pixels
onePlaneTolerance = 25;20;

radialDistToStim=zeros([size(stimCoM,1) numCells]);
axialDistToStim = zeros([size(stimCoM,1) numCells]);
StimDistance = zeros([size(stimCoM,1) numCells]);
for i=1:size(stimCoM,1);
    for k=1:numCells;
        D = sqrt(sum((stimCoM(i,:)-allCoM(k,:)).^2));
        radialDistToStim(i,k)=D;
        z = stimDepth(i)-allDepth(k);
        axialDistToStim(i,k) = z;
        StimDistance(i,k) = sqrt(sum((stimLoc(i,:)-allLoc(k,:)).^2));
        
    end
end

offTargetRisk = zeros([numel(roisTargets) numCells]);
for i=1:numel(roisTargets)
    Tg = roisTargets{i};
    
    if numel(Tg) == 1
        temp = radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0;
        temp2 = radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1;
    else
        temp = any(radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0);
        temp2 = any(radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1);
    end
    offTargetRisk(i,:) = temp | temp2;
end

%% Response to Vis stim

data = zdfData  ; %allSPData;    DFFdata; dfData; zData;

cellsToPlot = ~any(offTargetRisk)  & ~ROIinArtifact' & ~ismember(cellList,[HoloTargets{:}]) & pVisT<0.01;% & allDepth'>2.5 ;%& pVis2<0.01; %&  1:numCells;% ismember(cellList,TargetCells) & ~ROIinArtifact';

% cellsToPlot =  ~ROIinArtifact' & ismember(cellList,[HoloTargets{:}]);%1:numCells;% ismember(cellList,TargetCells) & ~ROIinArtifact';


trialsToPlot = lowMotionTrials & ~highRunTrials;

if ischar(trialsToPlot)
    trialsToPlot = ones(size(stimID));
end

disp([num2str(sum(cellsToPlot)) ' Cells. in ' num2str(sum(trialsToPlot)) ' Trials'])

baseline =1;
% normalize = 1;
% lim=([-0.5 0.75]);
lim=([-0.1 1]);

% trialsToPlot =

plotLog =[];
c=0;
X=[];
clear ax
figure(304);clf
for i = 1:numel(uniqueVisStim)
    v = uniqueVisStim(i);
    for k = 1:numel(uniqueStims)
        c=c+1;
        s=uniqueStims(k);
        
        ax(k,i) = subplot(2*numel(uniqueStims),numel(uniqueVisStim),i+2*numel(uniqueVisStim)*(k-1));
        x = data(cellsToPlot,:,visID==v & stimID == s & trialsToPlot);
        X{i,k}= x;
        
        m = mean(x,3);
        if baseline
         b = mean(m(:,1:3),2);
         m = m-b;
        end

        fillPlot(m,[],'ci');
        line([strt/1000*FR strt/1000*FR],get(gca,'ylim'));
        line([(strt+dura)/1000*FR (strt+dura)/1000*FR],get(gca,'ylim'));
        ylim(lim/3);
        xlim([0 30])
        
        
        xticks([0:FR:30])
        if k<numel(uniqueStims)
           xticklabels({''})
        else
            xticklabels({''})
%             xticklabels({'0' '1' '2' '3' '4' '5'})
        end
        
        if i~=1
            yticklabels({})
        end
%         
plotLog(c,:) = [c s v];
%         visR(c,:) = mean(m(:,win(1):win(2)),2);
        
%         if i==1 && k==1
%             ylabel('No Stim')
%             title('No Vis')
%         elseif i==1 && k==2
%             ylabel('Stim')
%         elseif i==2 && k==1
%             title('Vis')
%         end

if i==1 && k==1
    ylabel('No Stim')
    title('No Vis')
elseif i==1
    ylabel(['Stim: ' num2str(k)]);
elseif k==1
    title(['Vis: ' num2str(i)]);
end
% title(num2str(v));
%         
        
        subplot(2*numel(uniqueStims),numel(uniqueVisStim),i+numel(uniqueVisStim)+2*numel(uniqueVisStim)*(k-1));
        imagesc(mean(x,3));
        caxis(lim);
                xticks([0:FR:30])

        if k<numel(uniqueStims)
           xticklabels({''})
        else
            xticklabels({'0' '1' '2' '3' '4' '5'})
            xlabel('Seconds')
        end
         if i~=1
            yticklabels({})
        end
%         colorbar;
    end
    linkaxes(ax);
    
end
colormap plasma

%%co plot
% 
% colorList = {rgb('gray') rgb('FireBrick') rgb('ForestGreen')  rgb('DarkMagenta')};
% colorList = {rgb('lightGray') rgb('pink') rgb('lightGreen')  rgb('lightBlue')... };
%    rgb('darkgray') rgb('salmon') rgb('limeGreen') rgb('DarkTurquoise')...
%     rgb('dimGray') rgb('firebrick') rgb('ForestGreen') rgb('DarkCyan') };

% 
% colorList = {rgb('lightGray') rgb('pink')   rgb('lightBlue')... };
%    rgb('darkgray') rgb('salmon')  rgb('DarkTurquoise')...
%     rgb('dimGray') rgb('firebrick')  rgb('DarkCyan') };

cl=[];
colorList = {rgb('gray')};
cl = colormap('plasma');
nColors = c;
ncs = round(linspace(1,size(cl,1),nColors));
for i=1:nColors;
    colorList{end+1} =  cl(ncs(i),:);
end



cIdx=[];
figure(305);clf;
c=0;
for i = 1:numel(uniqueVisStim)
    v = visID(i);
    for k = 1:numel(uniqueStims)
        s=uniqueStims(k);
        c=c+1;
        cIdx(i,k) = c;
        
        x = X{i,k};
        m = mean(x,3);
        if baseline
        b = mean(m(:,1:3)');
        m = m-b';
        end

        fillPlot(m,[],'se',colorList{c},'none',colorList{c},0.5);
        
        line([strt/1000*FR strt/1000*FR],get(gca,'ylim'))
        line([(strt+dura)/1000*FR (strt+dura)/1000*FR],get(gca,'ylim'))
        
        line([visStart*FR visStart*FR],get(gca,'ylim'))
        line([visStop*FR visStop*FR],get(gca,'ylim'))
        
        xticks([0:6:30])
        xticklabels({'0' '1' '2' '3' '4' '5'})
        
        xlabel('Seconds')
        ylabel('Mean Z-Scored Fluorescence')
    end
end
title('All Conditions Superimposed')

% 
% visPlotList = {[1 2 3] [1 2 3] [1 2 3] [1] [1 2 3]};
% stimPlotList = {[1 2] [1 3] [1 4] [1 2 4] [1]};

% 
% visPlotList = {[1 2 3] [1 2 3]  [1] [1 2 3]};
% stimPlotList = {[1 2] [1 3]  [1 2 3] [1]};
%% ColorList Maker
cl=[];
colorList = {rgb('gray')};
cl = colormap('plasma');
nColors = numel(unique(visID));
ncs = round(linspace(1,size(cl,1),nColors));
for i=1:nColors;
    if i==3 | i==6
        colorList{end+1} = rgb('teal');
    else
        colorList{end+1} =  rgb('firebrick');
    end
end
colorList = {rgb('gray'), rgb('firebrick'), rgb('teal'), rgb('firebrick'),  rgb('teal'), rgb('firebrick'), rgb('teal'), rgb('purple')};

% nAlphas = numel(uniqueStims);
 alphaList = [0.75 0.5 0.25 0.75 0.5 0.25 0.75 0.5 0.25];%linspace(0.1,1,nAlphas);
 
 cl=[];
colorList2 = {rgb('gray')};
cl = colormap('viridis');
nColors = numel(unique(stimID));
ncs = round(linspace(1,size(cl,1),nColors));
for i=1:nColors;
    if i==3 | i==6
        colorList2{end+1} = rgb('teal');
    else
        colorList2{end+1} =  rgb('firebrick');
    end
end
colorList2 = {rgb('gray'), rgb('firebrick'), rgb('teal'), rgb('firebrick'),  rgb('teal'), rgb('firebrick'), rgb('teal'), rgb('plum')};

%% Fancy Plot Maker
stimNames = {'none', '3 cells tuned for 270', 'random 3 cells', '20 cells tuned for 270', 'random 20 cells',};
stimNames = {'Tuned Cells', 'Random Cells', '12 cell Ensembles', '6 Cell Ensembles',};

% visNames= {'none', '180 low contrast', '270 low contrast', 'none', '180 high contrast', '270 high contrast'}
visNames= {'none', '90\circ 50% Contrast', '135\circ 50% Contrast'}
visNames= {'none', '0\circ' '45\circ' '90\circ' '135\circ' '180\circ' '225\circ' '270\circ' '315\circ'}

%visPlotList =  {[2] [2] [3] [3]};% [2] [3] };%{[1 2 3] [ 1 2 3] [1] [2] [3]};
%stimPlotList = {[1 2 5] [1 3 6] [1 2 5] [1 3 6]};% [2] [1 6 7 8 ] [1 6 7 8] [1 6 7 8]};

%visPlotList = {[6] [6] [6] [6]};
%stimPlotList = {[4 1 8 2] [4 1 8 3] [7 1 8 5] [7 1 8 6]};
% stimPlotList = {[1 3 2] [1 3 2] [1 5 4] [1 5 4]};
stimPlotList = {[1 2 3] [1 4 5] [1 2 4] [1 3 5]};
stimPlotList ={[1 2 ]}; %L4 to L23
% stimPlotList = {[1 3 ]}; %L23 to L23
[ro co] = splitPretty(numel(stimPlotList),1,1);
visToTabThrough = 2; 8;1:numel(unique(visID));
% visToTabThrough = 8; %L23 to L23

for i=visToTabThrough;1:numel(unique(visID));
    visPlotList = repmat(i,numel(stimPlotList),1);
    figure(41);clf;
    for F = 1:numel(visPlotList);
        hold on
        subplot(ro,co,F);
        visToPlot = visPlotList(F);%[1 2 3];
        stimToPlot = stimPlotList{F};%[1 4];
        
        
        for L = 1:numel(visToPlot)
            v = visToPlot(L);
            for k = 1:numel(stimToPlot)
                si = stimToPlot(k);
                s=uniqueStims(si);
                c = cIdx(v,si);
                
                x = X{v,si};
                m = mean(x,3);
                if baseline
                    b = mean(m(:,1:3)');
                    m = m-b';
                end
                
                %         fillPlot(m,[],'se',colorList{c},'none',colorList{c},0.5);
                %         fillPlot(m,[],'se',colorList{v},'none',colorList{v},alphaList(si));
                fillPlot(m,[],'ci',colorList2{si},'none',colorList{si},0.25);
%                 line([strt/1000*FR strt/1000*FR],get(gca,'ylim'));
%                 line([(strt+dura)/1000*FR (strt+dura)/1000*FR],get(gca,'ylim'));
                yRange = get(gca,'ylim'); 
%                 line([visStart*FR visStart*FR],get(gca,'ylim'));
%                 line([visStop*FR visStop*FR],get(gca,'ylim'));
                   line([visStart*FR visStart*FR],yRange);
                line([visStop*FR visStop*FR],yRange);                             

                xticks([0:6:30]);
                xticklabels({'0' '1' '2' '3' '4' '5'});
                xlabel('Seconds');
                ylabel('Mean Z-Scored Fluorescence');
%                 title([stimNames(F) visNames(visToPlot)]);
            end
            
            
        end
        
    end
    if i<numel(visToTabThrough)
        
    pause
    end
end