out = All(2).out;
dfsource = 'zdfData';


%%%---- reliable visual responders ----%%%

vs = unique(out.vis.visID);


          
cellsToUse = ~out.anal.ROIinArtifact' &...
             out.anal.pVisR < 0.05;
         
for i=1:numel(vs)
    v=vs(i);
    
    trialsToUse = out.vis.visID == v &...
                  out.vis.lowMotionTrials;
              
              
    % split trials in half-ish
    ntrials = sum(trialsToUse);
    all_trials = find(trialsToUse);
    first_half = randsample(all_trials, floor(ntrials/2));
    second_half = all_trials(~ismember(all_trials, first_half));
    second_half = second_half(2:end);
    
    
    win = 12:18;
    dat1 = squeeze(mean(out.vis.(dfsource)(cellsToUse,win,first_half),2));
    bdat = squeeze(mean(out.vis.(dfsource)(cellsToUse,1:5,first_half),2));
    dat1 = dat1 - bdat;
    
    dat2 = squeeze(mean(out.vis.(dfsource)(cellsToUse,win,second_half),2));
    bdat = squeeze(mean(out.vis.(dfsource)(cellsToUse,1:5,second_half),2));
    dat2 = dat2 - bdat;
    
    % thsis gives a matrix..... I want pairwise
    c = corr(dat1',dat2');
    
    
    
end

         
         




% figure(71)
% clf
