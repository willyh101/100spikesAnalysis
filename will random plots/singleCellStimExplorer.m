clear redResps

figure(107)
clf
brws=5;
bcls=5;

subplotnum = 0;
for ind = 1:numel(All)
    out = All(ind).out;
    out.info
    
    
    dfsource = 'zdfData';
    
    if strcmp(dfsource, 'zdfData')
        dflbl = 'zdf';
    end
    % else
    %     dflbl = 'df/f';
    %     [dfData, ~] = computeDFFwithMovingBaseline(out.vis.allData);
    %     out.vis.dfData = dfData;
    %     [dfData, ~] = computeDFFwithMovingBaseline(out.exp.allData);
    %     out.exp.dfData = dfData;
    % end
    
    for stim = 2:numel(unique(out.exp.stimID))
        
        
        % ind  = 1;
        % stim = 4;
        vis  = [1]; % these should be the blank conditions
        
        
        
        us = unique(out.exp.stimID);
        s = us(stim);
        
        vs = unique(out.exp.visID);
        
        if max(vis) > max(vs)
            continue
        end
        
        v = vs(vis);
        
        h = out.exp.stimParams.roi{stim};
        tg = out.exp.holoTargets{h};
        tg(isnan(tg))=[];
        
        if numel(tg)~=1
            continue
        end
        
        
        
        
        
        
        %%------ stimmed cell figure --------%%
        % first, plot tuning curve, holo response, and vis resonse for the stimmed cell
        figure(66)
        clf
        sgtitle(['Stim ID ' num2str(stim) ', ROI ' num2str(tg)])
        
        % holo stim response
        subplot(1,3,1)
        
        trialsToUse = ismember(out.exp.stimID, s) &...
            ismember(out.exp.visID, v) &...
            out.exp.lowMotionTrials &...
            out.exp.stimSuccessTrial;
        
        if sum(trialsToUse) == 1
            continue
        end
        
        dat = squeeze(out.exp.(dfsource)(tg,:,trialsToUse));
        dat = dat-mean(dat(1:5,:),1);
        ntrials = sum(trialsToUse);
        
        [hL, hF] = fillPlot(dat',[],'ci');
        tcks = 0:6:30;
        xticks(tcks)
        xticklabels((tcks/6)-1)
        xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
        xlabel('Time')
        ylabel(dflbl)
        title('Holo stim response (40p, 40 Hz)')
        
        % vis response curve
        subplot(1,3,2)
        
        trialsToUse = out.vis.visID > 1 &...
            out.vis.lowMotionTrials;
        
        dat = squeeze(out.vis.(dfsource)(tg,:,trialsToUse));
        dat = dat-mean(dat(1:5,:),1);
        ntrials = size(dat,2);
        
        [hL, hF] = fillPlot(dat',[],'ci');
        tcks = 0:6:30;
        xticks(tcks)
        xticklabels((tcks/6)-1)
        xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
        xlabel('Time')
        ylabel(dflbl)
        title('Mean vis stim response, ori epoch')
        
        % tuning curve
        subplot(1,3,3)
        e = errorbar(out.anal.oriCurve(2:end,tg), out.anal.oriCurveSEM(2:end,tg));
        e.LineWidth = 2;
        e.Color = 'k';
        title(['Ori Tuning Curve, OSI=' num2str(out.anal.osi(tg))])
        ylabel(dflbl)
        xticks(1:12)
        xticklabels(0:30:330)
        xtickangle(45)
        xlabel('Direction')
        po = out.anal.prefOri(tg)-1;
        lbl_y = e.YData(po) + 1.1*e.YPositiveDelta(po);
        text(po, lbl_y, '*', 'Color','r','FontSize',30)
        
        
        
        
        %%------- other cells figure --------%%
        f67 = figure(67);
        clf
        rws = 1;
        cls = 5;
        
        sgtitle(['Stim ID ' num2str(stim) ', ROI ' num2str(tg)])
        
        
        
        % mean pop response all others
        subplot(rws,cls,1)
        
        trialsToUse = ismember(out.exp.stimID, s) &...
            ismember(out.exp.visID, v) &...
            out.exp.lowMotionTrials &...
            out.exp.stimSuccessTrial;
        
        cellsToUse = ~out.anal.ROIinArtifact' &...
            ~out.anal.offTargetRisk(stim-1,:) &...
                      out.anal.pVisR < 0.05;
        
        
        isred = out.red.isRed(cellsToUse);
        
        disp(['Total blue cells:     ' num2str(sum(out.red.isRed))])
        disp(['Blue cells vis resp:  ' num2str(sum(out.red.isRed & out.anal.pVisR < 0.05))])
        disp(['Blue cells used:      ' num2str(sum(isred))])
        
        dat = squeeze(mean(out.exp.(dfsource)(cellsToUse,:,trialsToUse),3));
        dat = dat-mean(dat(:,1:5),2);
        % popdat = squeeze(mean(dat,1));
        % ntrials = size(dat,2);
        
        % fillPlot is data,[],lineCol,edgeCol,faceCol,fAlpha
        
        [hL2, hF2] = fillPlot(dat(isred,:),[], 'ci', rgb('red'), [], rgb('red'), 0.5);
        hold on
        [hL1, hF1] = fillPlot(dat(~isred,:),[],'ci', rgb('black'), [], rgb('black'), 0.5);
        
        tcks = 0:6:30;
        xticks(tcks)
        xticklabels((tcks/6)-1)
        xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
        xlabel('Time')
        ylabel(dflbl)
        title('Mean Population response')
        legend([hL1,hL2], {'Pyramidal Cells','SST Cells'})
        
        % update big stim plot
        figure(107)
        subplotnum = subplotnum+1;
        subplot(brws, bcls, subplotnum)
        [hL2, hF2] = fillPlot(dat(isred,:),[], 'ci', rgb('blue'), [], rgb('blue'), 0.5);
        hold on
        [hL1, hF1] = fillPlot(dat(~isred,:),[],'ci', rgb('black'), [], rgb('black'), 0.5);
        
        tcks = 0:6:30;
        xticks(tcks)
        xticklabels((tcks/6)-1)
        xline(6,'--', {'Stim','Start'},'LineWidth',1,'LabelOrientation','horizontal','LabelHorizontalAlignment','center')
        xlabel('Time')
        ylabel(dflbl)
%         title('Mean Population response')
%         legend([hL1,hL2], {'Pyramidal Cells','SST Cells'})
        
        
        % scatter plot of all cells - no stim condition
        figure(67)
        subplot(rws,cls,2)
        
        trialsToUse = out.exp.stimID == us(1) &...
            ismember(out.exp.visID, [1 2]) &...
            out.exp.lowMotionTrials;
        
        dat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
        bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
        dat = dat - bdat;
        
        % dat = sort(dat);
        s1 = scatter(find(~isred),dat(~isred), 'filled');
        s1.MarkerFaceAlpha = 0.7;
        s1.MarkerFaceColor = rgb('grey');
        hold on
        s2 = scatter(find(isred),dat(isred), 'filled');
        s2.MarkerFaceColor = rgb('lightsalmon');
        
        yline(0, '--','LineWidth',1)
        ylabel(dflbl)
        % ax = gca();
        % ax.XAxis.Visible = 'off';
        xlabel('Cells')
        title({'Individual cell responses','No stim'})
        legend([s1,s2], {'Pyramidal Cells','SST Cells'})
        clear ax
        ax(1) = gca();
        
        
        % scatter plot of all cells
        subplot(rws,cls,3)
        
        trialsToUse = ismember(out.exp.stimID, s) &...
            ismember(out.exp.visID, v) &...
            out.exp.lowMotionTrials &...
            out.exp.stimSuccessTrial;
        
        dat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
        bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
        dat = dat - bdat;
        
        % dat = sort(dat);
        s1 = scatter(find(~isred),dat(~isred), 'filled');
        s1.MarkerFaceAlpha = 0.7;
        % s1.MarkerFaceColor = rgb('grey');
        hold on
        s2 = scatter(find(isred),dat(isred), 'filled');
        s2.MarkerFaceColor = 'r';
        
        yline(0, '--','LineWidth',1)
        ylabel(dflbl)
        % ax = gca();
        % ax.XAxis.Visible = 'off';
        xlabel('Cells')
        title('Individual cell responses')
        legend([s1,s2], {'Pyramidal Cells','SST Cells'})
        ax(2) = gca();
        linkaxes(ax)
        
        
        
        % % indiv cells by tuning (errorbar)
        % subplot(rws,cls,3)
        % colormap(f67, 'viridis')
        % oris = 0:30:330;
        % pref = oris(po);
        % all_pref = arrayfun(@(x) oris(x), out.anal.prefOri(out.anal.prefOri>1)-1);
        % diffToPO = mod(all_pref - pref,360);
        % dat2 = mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2);
        % error
        %
        %
        % s1 = scatter(1:sum(isred),sort(dat(isred)),out.red.redOSI*50, diffToPO, 'filled');
        % c = colorbar;
        % c.Label.String = 'PO difference';
        % c.Ticks = 0:30:330;
        % ylabel(dflbl)
        % xlabel('SST cells (sorted)')
        % title('SST response by difference to stim PO')
        
        
        % scatter plot difference
        
        subplot(rws,cls,4)
        
        
        trialsToUse = out.exp.stimID == us(1) &...
            ismember(out.exp.visID, v) &...
            out.exp.lowMotionTrials;
        
        dat1 = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
        bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
        dat1 = dat1 - bdat;
        
        trialsToUse = ismember(out.exp.stimID, s) &...
            ismember(out.exp.visID, v) &...
            out.exp.lowMotionTrials &...
            out.exp.stimSuccessTrial;
        
        dat2 = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
        bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
        dat2 = dat2 - bdat;
        
        new_dat = dat2 - dat1;
        
        s1 = scatter(find(~isred),new_dat(~isred), 'filled');
        s1.MarkerFaceAlpha = 0.5;
        % s1.MarkerFaceColor = rgb('blue');
        hold on
        s2 = scatter(find(isred),new_dat(isred), 'filled');
        s2.MarkerFaceColor = rgb('red');
        
        yline(0, '--','LineWidth',1)
        ylabel(dflbl)
        % ax = gca();
        % ax.XAxis.Visible = 'off';
        xlabel('Cells')
        title('Difference')
        legend([s1, s2], {'Pyramidal Cells', 'SST Cells'})
        ax(3) = gca();
        linkaxes(ax)
        
        
        % indiv cells by tuning (scatter since there aren't many)
        subplot(rws,cls,5)
        colormap(f67, 'viridis')
        oris = [nan 0:30:330];
        pref = oris(po);
        
        redOSI = All(ind).out.anal.osi(isred);
        redPO = All(ind).out.anal.prefOri(isred);
        redPO = arrayfun(@(x) oris(x), redPO);
        diffToPO = mod(redPO - pref,360);
        
        s1 = scatter(1:sum(isred),sort(dat(isred))',redOSI*100, diffToPO, 'filled');
        c = colorbar;
        c.Label.String = 'PO difference';
        c.Ticks = 0:30:330;
        ylabel(dflbl)
        xlabel('SST cells (sorted)')
        title('SST response by difference to stim PO')
        
        
        %%----- single cell scatter by sorted distance -----%%
        figure(87)
        clf
        
        rws = 2;
        cls = 2;
        
        % pyramids
        subplot(2,1,1)
        
        cellsToUse = ~out.anal.ROIinArtifact' &...
            ~out.anal.offTargetRisk(stim-1,:) &...
                    out.anal.pVisR < 0.05;
        
        trialsToUse = ismember(out.exp.stimID, us(1)) &...
            ismember(out.exp.visID, v) &...
            out.exp.lowMotionTrials &...
            out.exp.stimSuccessTrial;
        
        cell_dists = All(ind).out.anal.minDistbyHolo(stim,cellsToUse);
        [sorted_cell_dists, sort_idx] = sort(cell_dists);
        
        dat1 = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
        bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
        dat1 = dat1 - bdat;
        
        dat3 = squeeze(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2));
        dat3 = dat3 - bdat;
        
        
        trialsToUse = ismember(out.exp.stimID, s) &...
            ismember(out.exp.visID, v) &...
            out.exp.lowMotionTrials &...
            out.exp.stimSuccessTrial;
        
        dat2 = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2),3));
        bdat = squeeze(mean(mean(out.exp.(dfsource)(cellsToUse,1:5,trialsToUse),2),3));
        dat2 = dat2 - bdat;
        
        dat4 = squeeze(mean(out.exp.(dfsource)(cellsToUse,12:18,trialsToUse),2));
        dat4 = dat4 - bdat;
        
        new_dat = dat2 - dat1;
        sorted_resps = new_dat(sort_idx);
        
        diff_dat = dat4 - dat1;
        
        ps = zeros(size(new_dat));
        hs = zeros(size(new_dat));
        for c=1:size(diff_dat,1)
            %     [hval, pval] = ttest(dat4(c,:), mean(dat3(c,:)));
            [hval, pval] = ttest2(dat4(c,:), dat3(c,:));
            %     [hval, pval] = ttest(diff_dat(c,:));
            ps(c) = pval;
            hs(c) = hval;
        end
        
        ps_sorted = ps(sort_idx);
        isred = out.red.isRed(cellsToUse);
        ps_sorted(isnan(ps_sorted)) = 1;
        
        hs = hs(sort_idx);
        hs(isnan(hs)) = 0;
        % pyramids
        % subplot(rws,cls,1)
        
        colormap(flipud(colormap('viridis')))
        % s1 = scatter(sorted_cell_dists(~isred(sort_idx)), sorted_resps(~isred(sort_idx)), [], ps_sorted(~isred(sort_idx))', 'filled');
        
        s1 = scatter(sorted_cell_dists(~isred(sort_idx) & hs'), sorted_resps(~isred(sort_idx) & hs'), [], ps_sorted(~isred(sort_idx) & hs')', 'filled');
        hold on
        s2 = scatter(sorted_cell_dists(~isred(sort_idx) & ~hs'), sorted_resps(~isred(sort_idx) & ~hs'), [], ps_sorted(~isred(sort_idx) & ~hs')', 'filled');
        s1.MarkerEdgeColor = 'k';
        
        colorbar
        caxis([0.05 1])
        % s1.MarkerFaceAlpha = 0.5;
        
        yline(0, '--','LineWidth',1)
        ylabel(dflbl)
        % ax = gca();
        % ax.XAxis.Visible = 'off';
        xlabel('Distance (\mum)')
        title('Difference')
        legend([s1], {'Pyramidal Cells'})
        
        
        % blues
        subplot(2,1,2)
        
        colormap(flipud(colormap('viridis')))
        % s1 = scatter(sorted_cell_dists(~isred(sort_idx)), sorted_resps(~isred(sort_idx)), [], ps_sorted(~isred(sort_idx))', 'filled');
        
        s1 = scatter(sorted_cell_dists(isred(sort_idx) & hs'), sorted_resps(isred(sort_idx) & hs'), [], ps_sorted(isred(sort_idx) & hs')', 'filled');
        hold on
        s2 = scatter(sorted_cell_dists(isred(sort_idx) & ~hs'), sorted_resps(isred(sort_idx) & ~hs'), [], ps_sorted(isred(sort_idx) & ~hs')', 'filled');
        s1.MarkerEdgeColor = 'k';
        
        colorbar
        caxis([0.05 1])
        % s1.MarkerFaceAlpha = 0.5;
        
        yline(0, '--','LineWidth',1)
        ylabel(dflbl)
        % ax = gca();
        % ax.XAxis.Visible = 'off';
        xlabel('Distance (\mum)')
        title('Difference')
        legend([s1], {'SST Cells'})
        
        pause
    end
end

%% mean response of visually responsive SSTs and pyramids

figure(97)
clf








