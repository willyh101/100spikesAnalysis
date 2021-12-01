function individual_cells_plot(dist_all,cellResp_all,opts,xlabel_name)

subplot(2,2,1); hold on;
plot(dist_all,cellResp_all,'.')
plot([min(dist_all) max(dist_all)],0*[min(dist_all) max(dist_all)],'k--')
set(gca,'fontsize',16)
xlabel(xlabel_name)
ylabel('Pop Response (Mean \DeltaF/F)')

subplot(2,2,2); hold on;
[y_avg, y_std, y_sem] = calAvg(dist_all, cellResp_all, opts.distBins, diff(opts.distBins(1:2)));

plot(opts.distBins, y_avg, 'r-', 'LineWidth',2)
plot([0 opts.distBins], [0 opts.distBins*0], 'k--', 'LineWidth',2)
patch([opts.distBins fliplr(opts.distBins)],...
    [y_avg-y_sem fliplr(y_avg+y_sem)], 'r','FaceAlpha',.3)
set(gca,'fontsize',16)
xlabel(xlabel_name)
ylabel('Binned Pop Response')


subplot(2,2,3); hold on;
plot(dist_all,cellResp_all,'.')
plot([min(dist_all) max(dist_all)],0*[min(dist_all) max(dist_all)],'k--')
set(gca,'fontsize',16)
xlabel(xlabel_name)
ylabel('Pop Response (Mean \DeltaF/F)')
if strcmp(opts.distType,'min')
    xlim([15 50])
elseif strcmp(opts.distType,'mean')
    xlim([50 150])
elseif strcmp(opts.distType,'conn dist')
    xlim([4 7])
elseif strcmp(opts.distType,'gauss')
    xlim([6 8.5])    
end

subplot(2,2,4); hold on;
[y_avg, y_std, y_sem] = calAvg(dist_all, cellResp_all, opts.distBins, diff(opts.distBins(1:2)));

plot(opts.distBins, y_avg, 'r-', 'LineWidth',2)
plot([0 opts.distBins], [0 opts.distBins*0], 'k--', 'LineWidth',2)
patch([opts.distBins fliplr(opts.distBins)],...
    [y_avg-y_sem fliplr(y_avg+y_sem)], 'r','FaceAlpha',.3)
set(gca,'fontsize',16)
xlabel(xlabel_name)
ylabel('Binned Pop Response')

if strcmp(opts.distType,'min')
    xlim([15 50])
elseif strcmp(opts.distType,'mean')
    xlim([50 150])
elseif strcmp(opts.distType,'conn dist')
    xlim([4 7])
elseif strcmp(opts.distType,'gauss')
    xlim([6 8.5])
end


end


function [y_avg, y_std, y_sem] = calAvg(x, y, x_avg, dx, num_th)
% if nargin == 4; num_th = 20; end
if nargin == 4; num_th = 1; end

    y_avg = nan*zeros(size(x_avg));
    y_std = nan*zeros(size(x_avg));
    y_sem = nan*zeros(size(x_avg));
    
    for i = 1:length(x_avg)
        xids = logical((x>(x_avg(i)-dx/2)) .* (x<(x_avg(i)+dx/2)));
        if sum(xids(:)) >= num_th
            y_avg(i) = nanmean(y(xids));
            y_std(i) = nanstd(y(xids));
            y_sem(i) = nanstd(y(xids)) / sqrt(length(xids));
        end
    end
end

