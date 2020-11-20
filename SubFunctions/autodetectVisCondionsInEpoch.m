function [All] = autodetectVisCondionsInEpoch(All)

for ind = 1:numel(All)
    % Determine Vis epoch mean responses
    baseline =0;
    
    visDat=[];
    if baseline
        visDat = All(ind).out.vis.rdata - All(ind).out.vis.bdata;
    else
        visDat = All(ind).out.vis.rdata;
    end
    
    
    trialsToUse = All(ind).out.vis.lowMotionTrials & All(ind).out.vis.lowRunTrials;
    
    vs = unique(All(ind).out.vis.visID);
    
    mVisResp=[];
    for i = 1:numel(vs)
        v = vs(i);
        mVisResp(:,i) = nanmean(visDat(:,trialsToUse & All(ind).out.vis.visID==v),2);
    end
    
    
    % Determine Exp vis Responses
    
    expDat=[];
    if baseline
        expDat = All(ind).out.exp.rdData - All(ind).out.exp.bdata;
    else
        expDat = All(ind).out.exp.rdData;
    end
    
    trialsToUse = All(ind).out.exp.lowMotionTrials & All(ind).out.exp.lowRunTrials;
    
    
    vs = unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    us = unique(All(ind).out.exp.stimID);
    
    
    mExpResp=[];
    for i = 1:numel(vs)
        v = vs(i);
        s = us(1);
        mExpResp(:,i) = nanmean(expDat(:,...
            trialsToUse & All(ind).out.exp.visID==v & All(ind).out.exp.stimID == s),2);
    end
    
    
    %%Compare
    
    % % [coeff,score,latent,tsquared,explained,mu] = pca(mVisResp');
    % [coeff,score,latent,tsquared,explained,mu] = pca(visDat');
    %
    %
    % figure(4);clf
    % temp1 = coeff' * mVisResp;
    % scatter3(temp1(1,:),temp1(2,:),temp1(3,:))
    % % scatter3(score(:,1),score(:,2),score(:,3))
    % hold on
    %
    % temp = coeff' * mExpResp;
    % scatter3(temp(1,:),temp(2,:),temp(3,:))
    
    sz = size(mExpResp);
    if sz(2)==1
        visCode =1;
    else
        visCode=[];
        predictL1=[1];
        predictL2=[1];
        for i = 2:sz(2)
            try
            test = mExpResp(:,i);
            resid = abs(mVisResp-test);
            L1 = sum(resid);
            L2 = sqrt(sum(resid.^2));
            predictL1(i) = find(L1==min(L1),1);
            predictL2(i) = find(L2==min(L2),1);
            catch
                disp('Error')
                predictL1(i)=1;
                predictL2(i)=1;
            end
        end
        visCode = predictL1;
    end
    visCode
    All(ind).out.anal.visCode = visCode;
end
% temp = arrayfun(@(x) x.out.anal.visCode,All(IndsUsed),'uniformoutput',0)

%% The Auto Detect approach worked but I checked it anyways
% All(9).out.anal.visCode = [1 4];
% All(12).out.anal.visCode = [1 5 7];
% All(15).out.anal.visCode = [1 6 9 6 9];
