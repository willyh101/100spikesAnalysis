function [outVars] = plotCorrelationResponse(All,outVars,stringCorrType)
numExps = numel(All);
%% Plot
    disp('Calculating...')

clear popResponseCorr
for ind = 1:numExps
    
    try
    corrToUse  = eval(['All(ind).out.anal.' stringCorrType]); %AllCorr;
    catch
        disp('That Corr type not available')
        return
    end
    
    vs =  unique(All(ind).out.exp.visID);
    vs(vs==0)=[];
    respMat = All(ind).out.anal.respMat;
    baseMat = All(ind).out.anal.baseMat;
    
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    
    clear popRespCorr minDistbyHolo cellsToUse
    for v = 1:numel(vs)
        for i= 1:numel(All(ind).out.exp.stimParams.Seq)
            holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
            if i==1
                cellsToUse = ~ROIinArtifact';
            else
                cellsToUse = ~ROIinArtifact'  & ~offTargetRisk(holo,:);
            end
            %             popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if i~=1
                Tg=All(ind).out.exp.holoTargets{holo};
                Tg(isnan(Tg))=[];
                
                distCorr = corrToUse(Tg,:);
                if numel(Tg)==1
                    minDist = distCorr;
                else
                    minDist = nanmean(distCorr);
                end
                
                if numel(Tg)==0
                    minDistbyHolo(i,:) = ones([1 size(minDist,2)])*1000;
                else
                    minDistbyHolo(i,:) = minDist;
                end
                distBins = linspace(-0.5,0.5,40);
                for d = 1:numel(distBins)-1
                    cellsToUse = ~ROIinArtifact' &...
                        ~offTargetRisk(holo,:) &...
                        minDist > distBins(d) &...
                        minDist <= distBins(d+1) ;
                    popRespCorr(i,v,d) = nanmean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                    
                    noHoloEquivalent = nanmean(squeeze(respMat(1,v,cellsToUse) - baseMat(1,v,cellsToUse)));
                    popRespCorrSub(i,v,d) =  popRespCorr(i,v,d) - noHoloEquivalent;
                end
            end
        end
    end
    popRespCorr(1,:,:)=[];
    popRespCorrSub(1,:,:)=[];
    
    popResponseCorr{ind} = popRespCorr;
    popResponseCorrSub{ind} = popRespCorrSub;
end


%
% popResponseCorr = cell2mat(popResponseCorr(:));
% popResponseCorr(numSpikesEachStim==0)=[];

outVars.popResponseCorr = popResponseCorr;
outVars.popResponseCorrSub = popResponseCorrSub;


temp = cellfun(@(x) permute(x(:,1,:),[1 3 2]),popResponseCorr,'uniformoutput',0) ;
% temp = cellfun(@(x) squeeze(x(:,end,:)),popResponseCorr,'uniformoutput',0) ;

% temp = cellfun(@(x) squeeze(x(:,end,:)),popResponseCorrSub,'uniformoutput',0) ;
% temp = cellfun(@(x) squeeze(x(:,1,:)),popResponseCorrSub,'uniformoutput',0) ;

% temp = cellfun(@(x) squeeze(x(:,round(end/2),:)),popResponseCorr,'uniformoutput',0) ;

EnsCorR = cat(1,temp{:});

figure(15);clf
% subplot(1,3,1)

ensemblesToUse = outVars.ensemblesToUse;
numCellsEachEns = outVars.numCellsEachEns;
uniqueEns = unique(numCellsEachEns(ensemblesToUse));

numEns = numel(uniqueEns);

colorList = colorMapPicker(numEns,'parula');

hold on
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse );
    data = EnsCorR(ens2plot,:);
    
    
    e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data)./sqrt(sum(~isnan(data))));
    
    %         e = errorbar(distBins(1:end-1),nanmean(data,1),nanstd(data));
    
    e.Color = colorList{i};
    e.LineWidth = 2;
    hack = num2cell(distBins(1:end-1));
    hack = cellfun(@(x) num2str(x),hack,'uniformoutput',0);
    %     p = plotSpread(data, 'xNames', hack, 'showMM', 4);
    % fancyPlotSpread(data,hack)
    %     names{i} = string(uniqueEns(i));
    %     avg(i) = mean(popResponseEns(ens2plot));
    %     err(i) = sem(popResponseEns(ens2plot));
    %     ns(i) = numel(popResponseEns(ens2plot));
end
xlim([-0.4 0.4])
legend(cellfun(@(x) num2str(x),num2cell(uniqueEns),'uniformoutput',0))
% legend('small', 'medium', 'large')

r = refline(0);
r.LineStyle = ':';
r.Color = rgb('grey');
r.LineWidth = 2;

ylabel('Pop Resp to HoloStim')
xlabel('Responder to Ensemble Correlation')