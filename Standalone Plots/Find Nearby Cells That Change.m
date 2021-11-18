
 c=1;cellPairOfInterest=[]; cellPairOfInterestEnsemble=[];
for ind = IndsUsed
    
    eligibleEns = find(outVars.ensemblesToUse & ensIndNumber==ind);
    
    HTs=[];
    for i = 1:numel(eligibleEns)
        e = eligibleEns(i);
        s = outVars.ensHNumber(e); %i think this is the order of resp
        h = All(ind).out.exp.stimParams.roi{s};
        htg = All(ind).out.exp.holoTargets{h};
        %        htg(isnan(htg))=[];
        HTs(i,:) = htg;
    end
    
    uHTs = unique(HTs);
    uHTs(isnan(uHTs))=[];
    
    countHTs = [];
    for i=1:numel(uHTs)
        countHTs(i) = sum(HTs(:)==uHTs(i));
    end
    
    targsHitMulti = uHTs(countHTs>1);
    
    allDist = All(ind).out.anal.allDistance;
    range =25; %distance away to consider
    adjacentCells =[];
    for i=1:numel(targsHitMulti)
        t = targsHitMulti(i);
        adjacentCells = find(allDist(t,:)<=range & allDist(t,:)>0);
        adjacentCellsAll{i} = adjacentCells;
        
        es = eligibleEns(any(HTs==t,2));
        
        if ~isempty(adjacentCells)
            for k=1:numel(adjacentCells)
                a = adjacentCells(k);
                AdjCellResp=[];
                for j = 1:numel(es)
                    e = es(j);
                    s = outVars.ensHNumber(e);
                    h = All(ind).out.exp.stimParams.roi{s};
                    htg = All(ind).out.exp.holoTargets{h};
                    if ismember(a,htg)
                        AdjCellResp(j) = NaN;
                    else
                        AdjCellResp(j)=squeeze(All(ind).out.anal.respMat(s,1,a)-All(ind).out.anal.baseMat(s,1,a));
                    end
                end
                %                                AdjCellResp
                upperBound = 0.55;
                if min(AdjCellResp)<0 && max(AdjCellResp)>upperBound;
                    cellPairOfInterest(c,:) = [ind t a];
                    cellPairOfInterestEnsemble{c} = [es(AdjCellResp<0) es(AdjCellResp>upperBound)];
                    c=c+1;
                end
            end
        end
    end
    
    
end

cellPairOfInterest

%%
figure(10001);clf
 for i =1:size(cellPairOfInterest,1)
     clf
     ind = cellPairOfInterest(i,1);
     t = cellPairOfInterest(i,2);
     a = cellPairOfInterest(i,3);
     
     e1 = cellPairOfInterestEnsemble{i}(1);
     e2 = cellPairOfInterestEnsemble{i}(2);

     try
     ax = subplot(2,6,[1 7]);
     [All(ind)] = plotCellMask(All(ind),t,ax);
     title(num2str(t))
     catch;end
     
     ax = subplot(2,6,2);
     plotCaRaster(All,outVars,e1,t,ax)
     ax = subplot(2,6,8);
     plotCaRaster(All,outVars,e2,t,ax)
     
     ax = subplot(2,6,3);
     plotCellMean(All,outVars,e1,t,rgb('darkmagenta'),ax)
     ax = subplot(2,6,9);
     plotCellMean(All,outVars,e2,t,rgb('darkmagenta'),ax)

     try
     ax = subplot(2,6,[4 10]);
     [All(ind)] = plotCellMask(All(ind),a,ax);
         title(num2str(a))
     catch; end
     
         ax = subplot(2,6,5);
     plotCaRaster(All,outVars,e1,a,ax)
     ax = subplot(2,6,11);
     plotCaRaster(All,outVars,e2,a,ax)

          ax = subplot(2,6,6);
     plotCellMean(All,outVars,e1,a,rgb('steelblue'),ax)
     ax = subplot(2,6,12);
     plotCellMean(All,outVars,e2,a,rgb('steelblue'),ax)
     
     pause
 end
 
 disp('done')
 
 %% Pull Example Images
 
 ind =31;
 cell1 = 390;
 cell2 = 386;
 depth = All(ind).out.exp.allDepth(390);
 nDepths = 3; 
 
 e1 =   822;
 e2  =  830;
 
 us = unique(All(ind).out.exp.stimID); 
 
 s1 = outVars.ensHNumber(e1); %i think this is the order of resp
 sID1 = All(ind).out.exp.uniqueStims(s1);
 h1 = All(ind).out.exp.stimParams.roi{s1};
 
 s2 = outVars.ensHNumber(e2); %i think this is the order of resp
 sID2 = All(ind).out.exp.uniqueStims(s2);
 h2 = All(ind).out.exp.stimParams.roi{s2};
 
%  trialsToUse1 = All(ind).out.exp.lowMotionTrials &...
%     All(ind).out.exp.lowRunTrials &...
%     All(ind).out.exp.stimSuccessTrial &...
%     All(ind).out.exp.stimID == us(s1) & ...
%     All(ind).out.exp.visID == 1;% | All(ind).out.exp.visID == 0);

 trials1 = All(ind).out.anal.defaultTrialsToUse & All(ind).out.exp.stimID==sID1;
 trials2 = All(ind).out.anal.defaultTrialsToUse & All(ind).out.exp.stimID==sID2;
 
 
 pathForFiles = 'T:\Ian\I154_2\210927\6';
 
 tic;
 d = dir(pathForFiles);
 d(1:2)=[];
 
 T = All(ind).out.exp.Tarray;
 numFrames = size(T{depth},1);
 
 disp('Ensemble 1')
 clear MCG
 trials1List = find(trials1);
numEl = numel(trials1List);
fprintf('|')
for i=1:numEl
 fprintf('-');
end
disp('|')

fprintf('|')
 for i=1:numEl;
     fprintf('.')
     
     t = trials1List(i);
     fN  = fullfile(d(t).folder,d(t).name); 
      v = ScanImageTiffReader(fN).data();
      v = permute(v,[2 1 3]);
      
      gg = v(:,:,1:2:end);
      gg = gg(:,:,depth:nDepths:end);

%       ss=size(gg);
%       mcg=zeros(ss);
%       for n=(1:min(ss(3),numFrames));
%           mcg(:,:,n)=circshift(gg(:,:,n),round(squeeze(T{depth}(n,t,:))));
%       end
      mcg = gg(:,:,1:numFrames); %simpleAlignTimeSeries(gg(:,:,1:numFrames),0);
      
      MCG(:,:,:,i) = mcg(:,:,1:numFrames);
 end
 
 sz  = size(MCG);
 MCUnr= reshape(MCG, [sz(1) sz(2) sz(3)*sz(4)]);
 [MCG dxs dys] =  simpleAlignTimeSeries(MCUnr,1);
 MCG1 = reshape(MCG, sz);

 fprintf('|')
 
 disp('Now Ensemble 2')
 
 clear MCG
 
  trials2List = find(trials2);
numEl = numel(trials2List);
fprintf('|')
for i=1:numEl
 fprintf('-');
end
disp('|')

fprintf('|')
 for i=1:numEl;
     fprintf('.')
     
     t = trials2List(i);
     fN  = fullfile(d(t).folder,d(t).name); 
      v = ScanImageTiffReader(fN).data();
      v = permute(v,[2 1 3]);
      
      gg = v(:,:,1:2:end);
      gg = gg(:,:,depth:nDepths:end);

%       ss=size(gg);
%       mcg=zeros(ss);
%       for n=(1:min(ss(3),numFrames));
%           mcg(:,:,n)=circshift(gg(:,:,n),round(squeeze(T{depth}(n,t,:))));
%       end
      mcg = gg(:,:,1:numFrames); %simpleAlignTimeSeries(gg(:,:,1:numFrames),0);
      
      MCG(:,:,:,i) = mcg(:,:,1:numFrames);
 end
 
 sz  = size(MCG);
 MCUnr= reshape(MCG, [sz(1) sz(2) sz(3)*sz(4)]);
 [MCG dxs2 dys2] =  simpleAlignTimeSeries(MCUnr,1);
 MCG2 = reshape(MCG, sz);
 
 
 fprintf('|')
disp('done')
 toc
 
 
 %%
 M1 = mean(MCG1,4);
 M2 = mean(MCG2,4);
 
 rcwin = All(ind).out.anal.recWinUsed;
 bwin = All(ind).out.anal.bwinToUse;
 
 M1R  = mean(M1(:,:,rcwin(1):rcwin(2)),3) ;%- mean(M1(:,:,bwin(1):bwin(2)),3);
 M2R  = mean(M2(:,:,rcwin(1):rcwin(2)),3) ;%- mean(M2(:,:,bwin(1):bwin(2)),3);

 figure(3);
 ax1 =  subplot(1,2,1)
%  imagesc(mean(mean(MCG1,4),3));
imagesc(M1R)
axis square
 ax2  = subplot(1,2,2);
%  imagesc(mean(mean(MCG2,4),3));
imagesc(M2R)
 linkaxes([ax1 ax2])
 axis square
 
  offsets = All(ind).out.exp.offsets; 

 xy1 = All(ind).out.exp.allCoM(cell1,:)-offsets;
 xy2 = All(ind).out.exp.allCoM(cell2,:)-offsets;
 
 
 mn = min([xy1; xy2]);
 mx = max([xy1; xy2]);
 
 xlim([mn(1)-10 mx(1)+10]);
 ylim([mn(2)-10 mx(2)+10]);
%  
subplot(1,2,1)
 circle(xy1(1)+mean(dxs),xy1(2)+mean(dys),5,rgb('magenta'));
  circle(xy2(1)+mean(dxs2),xy2(2)+mean(dys2),5,rgb('steelblue'))
  
  subplot(1,2,2)
   circle(xy1(2)+mean(dys),xy1(1)+mean(dxs),5,rgb('magenta'));
  circle(xy2(2)+mean(dys2),xy2(1)+mean(dxs2),5,rgb('steelblue'))
