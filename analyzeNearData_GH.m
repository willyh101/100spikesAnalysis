%%
% Designed to analyze the data corresponding cells closest to the ensemble 
% (though the data could be any distance away, one would just have to
% adjust the type of t-test done below)
%
% nearToEnsData: cell array containing the population response for those 
%                closest the the ensemble, for the different bins
% bins2: the bin edges for criteria 2 (e.g., 1- circ var PO)
%
% If bins2 == 4 (meaning three bins), it only plots the first and last row
%
% Performs a one-sided t-test for all pairs
% Hypothesis: As you go from closest cluster to far apart, you will see 
% an increase in pop. response
%
% Can be adjusted to be two-sided, and will test to see if the means are
% different (see ttest2 below)
%
% Written by Gregory Handy, (Unviersity of Chicago, 2021)
%%
function [] = analyzeNearData_GH(nearToEnsData,bins2)

% Plot the distributions for the first data point
figure(); clf

bins_to_test{1} = [1 2 3];
if numel(bins2) == 4
    bins_to_test{2} = [7 8 9];
    bins_to_test{3} = [7 1];
else
    bins_to_test{2} = [4 5 6];
    bins_to_test{3} = [4 1];
end

num_bins_to_test = 3;
for hh =1:num_bins_to_test

    ymax(hh) = 0;
    subplot(num_bins_to_test,1,hh)
    hold off
    for ii = 1:length(bins_to_test{hh})
        temp1 = nearToEnsData{bins_to_test{hh}(ii)};
        ymax(hh) = max(ymax(hh),max(temp1));
        plot(ii+0*temp1-0.05+0.1*rand(length(temp1),1),temp1,'.','markersize',16)
        hold on
        plot([ii-0.25:0.01:ii+0.25],mean(temp1)+0*[ii-0.25:0.01:ii+0.25],'k-','linewidth',1.5)
    end
    
    set(gca,'fontsize',16)
    xticks([1 2 3])
    xticklabels({'Close Together' 'Medium' 'Far Apart'})
    if hh == 1 
        title('Untuned')
    elseif hh == 2
        title('Co-tuned')
    else
        subplot(3,1,3)
        xticks([1 2])
        xticklabels({'Co-tuned', 'Untuned'})
        title('Close together')
    end
end

% Perform one-sided (or two-sided) t-tests and plot stars
ttest_pairs{1} = [1 2; 1 3; 2 3];

if numel(bins2) == 4
    ttest_pairs{2} = [7 8; 7 9; 8 9];
    ttest_pairs{3} = [7 1];
else
    ttest_pairs{2} = [4 5; 4 6; 5 6];
    ttest_pairs{3} = [4 1];
end


linesizes = [1 2; 1 3; 2 3];
for hh = 1:num_bins_to_test
    subplot(3,1,hh)
    
    inner_count = 1;
    for ii = 1:size(ttest_pairs{hh},1)
        [~, p_val] = ttest2(nearToEnsData{ttest_pairs{hh}(ii,1)},nearToEnsData{ttest_pairs{hh}(ii,2)},'tail','left');
        plot(linesizes(ii,:), [ymax(hh)+0.03+.2*(ii-1) ymax(hh)+0.03+.2*(ii-1)],'k','linewidth',2)
        
        txt = sprintf('p=%.3f',p_val);
        text(mean(linesizes(ii,:))-.3,ymax(hh)+0.1+.2*(ii-1),txt,'fontsize',16)
        if p_val < 0.05
            plot(mean(linesizes(ii,:)), ymax(hh)+0.1+.2*(ii-1),'k*','markersize',10)
        end
        
        inner_count = inner_count + 1;
    end
    
    ylim_temp = ylim;
    ylim([ylim_temp(1) ymax(hh)+0.1+.2*(ii-1)+0.1])
    
end


end

