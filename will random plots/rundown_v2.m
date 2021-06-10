figure(77); clf

out = All(1).out;
dfsource = 'zdfData';

vis = [1];

rws = numel(unique(out.exp.stimID))-1;

c = 0;
for i = 2:numel(unique(out.exp.stimID))
    
    c = c+1;
    
    us = unique(out.exp.stimID);
    s = us(i);

    vs = unique(out.exp.visID);
    v = vs(vis);
    
    h = out.exp.stimParams.roi{i};
    tg = out.exp.holoTargets{h};
    tg(isnan(tg))=[];
    
    trialsToUse = ismember(out.exp.stimID, s) &...
              ismember(out.exp.visID, v) &...
              out.exp.lowMotionTrials &...
              out.exp.stimSuccessTrial;
          
    if sum(trialsToUse) == 0
        disp(['skipped stim ' num2str(s)])
        c = c - 1;
        continue
    end

    dat = squeeze(out.exp.(dfsource)(tg,:,trialsToUse));
    
    if ndims(dat) > 2
        dat = squeeze(mean(dat, 1));
    end
    
    dat = dat-mean(dat(1:5,:),1);

    ntrials = sum(trialsToUse);

    resp = mean(dat(12:18,:),1);
    
    %%----- Trialwise response plot -----%%
        
    subplot(rws,2,c)
    
    scatter(find(trialsToUse), resp, 'filled', 'MarkerFaceAlpha', 0.7)

    hold on
    
    ft = polyfit(find(trialsToUse), resp, 1);
    ft_xs = 1:length(trialsToUse);
    ftvals = polyval(ft, ft_xs);
    plot(ft_xs, ftvals, 'LineWidth', 1, 'Color', 'k')
    title(['Cell ' num2str(i-1)])
    
    %%----- Trace response plot -----%%
    
    c = c+1;
    subplot(rws,2,c)
    
    [hL, hF] = fillPlot(dat',[],'ci');
    tcks = 0:6:30;
    xticks(tcks)
    xticklabels((tcks/6)-1)
    xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
    xlabel('Time')
    ylabel(dflbl)
    title('Holo stim response (40p, 40 Hz)')
    
end
