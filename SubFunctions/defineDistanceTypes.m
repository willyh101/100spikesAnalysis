function [All, outVars] = defineDistanceTypes(All, outVars)
numExps = numel(All);

for ind = 1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    if isfield(All(ind).out.exp, 'dataToUse')
        dataToUse = All(ind).out.exp.dataToUse;
    else
        disp(['ind ' num2str(ind) '. no data to use, using zdfData']);
        dataToUse = All(ind).out.exp.zdfData;
    end
    
    
    stimCoM = All(ind).out.exp.stimCoM;
    numCells = size(dataToUse,1);
    allCoM = All(ind).out.exp.allCoM;
    stimDepth = All(ind).out.exp.stimDepth;
    allDepth = All(ind).out.exp.allDepth;
    muPerPx = 800/512;
    
    allLoc = [allCoM*muPerPx (allDepth-1)*30];
    stimLoc = [stimCoM*muPerPx (stimDepth-1)*30];

    StimDistance = zeros([size(stimCoM,1) numCells]);
    for i=1:size(stimCoM,1);
        for k=1:numCells;
            StimDistance(i,k) = sqrt(sum((stimLoc(i,:)-allLoc(k,:)).^2));
        end
    end
    All(ind).out.anal.StimDistance = StimDistance;
    
    allDistance = zeros([numCells numCells]);
    for i=1:numCells;
        for k=1:numCells;
            allDistance(i,k) = sqrt(sum((allLoc(i,:)-allLoc(k,:)).^2));
        end
    end
    All(ind).out.anal.allDistance = allDistance;
    
    roiDistance = zeros([size(stimCoM,1) size(stimCoM,1)]);
    for i=1:size(stimCoM,1);
        for k=1:size(stimCoM,1);
            roiDistance(i,k) = sqrt(sum((stimLoc(i,:)-stimLoc(k,:)).^2));
        end
    end
    All(ind).out.anal.roiDistance = roiDistance;
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% Determine the Ensemble CoCorrelation


ensMeaD=[];ensGeoD=[];ensMaxD=[];ensMinD=[];
for ind = 1:numExps
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    
    clear ensembleMeanDist ensembleGeoMeanDist ensembleMaxDist ensembleMinDist
    %     for i =1:numel(All(ind).out.exp.holoTargets)
    %         ht = All(ind).out.exp.holoTargets{i};
    %         ht(isnan(ht))=[];
    c=0;
    for i= 1:numStims
        holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
        if holo>0
            c=c+1;
            ht = All(ind).out.exp.holoTargets{holo};
            ht(isnan(ht))=[];
            
            rt = All(ind).out.exp.rois{holo};
            

            distToUse = All(ind).out.anal.roiDistance;
            distMat = distToUse(rt,rt);
            distMat(logical(eye(numel(rt))))=nan;
            
            ensembleMeanDist(c)     = nanmean(distMat(:));
            tempDist = distMat(~isnan(distMat));
            try
                ensembleGeoMeanDist(c)  = geo_mean(tempDist(:));
            catch
                ensembleGeoMeanDist(c)  = geomean(tempDist(:));
            end
            ensembleMaxDist(c)      = max(distMat(:));
            ensembleMinDist(c)      = min(distMat(:));
            
        end
    end
    All(ind).out.anal.ensembleMeanDist      = ensembleMeanDist;
    All(ind).out.anal.ensembleGeoMeanDist   = ensembleGeoMeanDist;
    All(ind).out.anal.ensembleMaxDist       = ensembleMaxDist;
    All(ind).out.anal.ensembleMinDist       = ensembleMinDist;
    
    ensMeaD = cat(2,ensMeaD,ensembleMeanDist);
    ensGeoD = cat(2,ensGeoD,ensembleGeoMeanDist);
    ensMaxD = cat(2,ensMaxD,ensembleMaxDist);
    ensMinD = cat(2,ensMinD,ensembleMinDist);
end

outVars.ensMeaD=ensMeaD;
outVars.ensGeoD=ensGeoD;
outVars.ensMaxD=ensMaxD;
outVars.ensMinD=ensMinD;


