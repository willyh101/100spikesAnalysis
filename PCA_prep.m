

%%

opts.runThreshold = 6 ; %trials with runspeed below this will be excluded
opts.runValPercent = 0.75; %percent of frames that need to be below run threshold
%%
out.vis.lowRunTrials = ones(1,size(All(1).out.vis.runVal,1));
for trial_num = 1:size(All(1).out.vis.runVal,1)
    
    percent_run = length(find(out.vis.runVal(trial_num,1:winToUse(2))<opts.runThreshold))...
        /length(out.vis.runVal(trial_num,1:winToUse(2)));
    
    if percent_run<=opts.runValPercent
        out.vis.lowRunTrials(trial_num) = 0;
    end
end

%%



% timeframes for the visual stimulus 
recWinSec=recWinRange+out.vis.visStart;
winToUse = min(round(recWinSec*out.info.FR),[inf sz(2)]);

bwinToUse = max(round([0 out.vis.visStart]*out.info.FR),[1 1]);

% neuron_id = ceil(size(All(1).out.vis.zdfData,1)*rand())
% find(All(1).out.anal.osi>0.8)
neuron_id = 14;

%visID: 2 is 0 degrees, and 4 is 90 degrees
visIDToTest = [2 4];

figure(1); clf;
figure(2); clf;
for ii = 1:length(visIDToTest)
    
    figure(1);
    test_indices = find(out.vis.visID==visIDToTest(ii) &...
        out.vis.lowMotionTrials &...
        out.vis.lowRunTrials);
    subplot(1,2,ii);hold on;
    temp = squeeze(out.vis.zdfData(neuron_id,:,test_indices));
    plot(temp)
    plot(mean(squeeze(out.vis.zdfData(neuron_id,:,test_indices)),2),'k','linewidth',1.5)
    
    plot(winToUse(1)+0*[min(min(temp)) max(max(temp))],[min(min(temp)) max(max(temp))],'k--')
    plot(winToUse(2)+0*[min(min(temp)) max(max(temp))],[min(min(temp)) max(max(temp))],'k--')
    
    
    ave_save(ii) = mean(mean(temp(winToUse(1):winToUse(2),:)));
    
    figure(2); hold on
    plot(mean(squeeze(out.vis.zdfData(neuron_id,:,test_indices)),2),'k','linewidth',1.5);
    plot(winToUse(1)+0*[min(min(temp)):max(max(temp))],[min(min(temp)):max(max(temp))],'k--')
    plot(winToUse(2)+0*[min(min(temp)):max(max(temp))],[min(min(temp)):max(max(temp))],'k--')
    
end


figure(3); clf; hold on;
plot(All(1).out.anal.oriCurve([2 4],neuron_id))
plot(ave_save)

%%

trialsToUse = All(1).out.vis.visID==2 &...
            All(1).out.vis.lowMotionTrials &...
            All(1).out.vis.lowRunTrials;
  
        
find(All(1).out.vis.lowRunTrials==0)


%%        
        
        
sum(trialsToUse)
test_indices = find(All(1).out.vis.visID==2);
length(test_indices)

%%
%%

% X = training_centered'*training_centered;
% [V_GH,D] = eig(X);
% 
% [B,I] = sort(diag(D),'descend');
% 
% D = diag(B);
% V_GH = V_GH(:,I);
% 
% D2(end)=0;
% sigma = zeros(62,634);
% sigma(1:62,1:62) = sqrt(D(1:62,1:62));
% 
% sigma(1,1)^2/(62-1)
% 
% plot(cumsum(diag(sigma).^2)/sum(diag(sigma).^2))
% 
% for ii = 1:size(training_centered,1)
%     U_GH_v2(:,ii) = training_centered*V_GH(:,ii)/sqrt(D(ii,ii));
% end
% 
% norm(training_centered-U_GH_v2*sigma*V_GH')
% 
% test_S = U_GH_v2*sigma;
% test_S = test_S(1:62,1:62);
% 
% test_C = V_GH(:,1:62);
% 
% 
% norm(training_centered-test_S*test_C')
% 
% %%
% [U,S,V] = svd(training_centered);
% 
% norm(training_centered-U*S*V')
% 
% %%
% S(1,1)
% LATENT(1,1)
% %%
% figure(3); clf; hold on;
% plot(All(1).out.anal.oriCurve([2 4],neuron_id))
% plot(ave_save)
% 
% %%
% 
% trialsToUse = All(1).out.vis.visID==2 &...
%             All(1).out.vis.lowMotionTrials &...
%             All(1).out.vis.lowRunTrials;
%   
%         
% find(All(1).out.vis.lowRunTrials==0)
% 
% 
% %%        
%         
%         
% sum(trialsToUse)
% test_indices = find(All(1).out.vis.visID==2);
% length(test_indices)
