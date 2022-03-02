muPerPix = 800/512;

disp('Rematching Targets')
for ind=1:numExps
   fprintf('.')
    
    stimDepth = All(ind).out.exp.stimDepth;
    allDepth = All(ind).out.exp.allDepth;
    
    stimCoM = All(ind).out.exp.stimCoM;
    allCoM = All(ind).out.exp.allCoM;
    
     numRois = numel(stimDepth);
    numCells = size(allCoM,1);
    
   radialStimDistance = nan([numRois numCells]); %in um
   for i=1:numRois;
      for k=1:numCells;
          if stimDepth(i)==allDepth(k);
             radialStimDistance(i,k) =  sqrt(sum((stimCoM(i,:) - allCoM(k,:)).^2))*muPerPix;
          end
      end
   end
    
[distToClosest newTargetedCells] = min(radialStimDistance,[],2);
   
newTargetedCells(distToClosest>opts.matchDistanceMicrons)=nan;
All(ind).out.exp.targetedCells = newTargetedCells;


for i =1:numel(All(ind).out.exp.rois);
    rt = All(ind).out.exp.rois{i};
    ht = newTargetedCells(rt);
%     ht(isnan(ht))=[];
    All(ind).out.exp.holoTargets{i}=ht';
end
   
end
disp('matched')