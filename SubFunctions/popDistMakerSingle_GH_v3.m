function [popRespDistEns, distVals,cellResp] = popDistMakerSingle_GH_v3(opts,All,CellToUseVar,plotFlag,ensID)
%assumes one All
%% make distance plots
distType = opts.distType;
distBins = opts.distBins; %[0:25:1000];
numDist =numel(distBins)-1;
numExps = numel(All);

c=0;
ind = 1;
%add additional constraints to cellsToUse
if isempty(CellToUseVar)
    cellToUseLimit = ones([1 size(All(ind).out.anal.respMat,3)]);
elseif islogical(CellToUseVar)
    cellToUseLimit=CellToUseVar;
else
    try
        cellToUseLimit = eval(['All(ind).out.' CellToUseVar]);
    catch
        disp('ERROR variable not found')
        return
    end
end
%First create stimDistance
% if ~isfield(All(ind).out.anal,'StimDistance') || strcmp(distType,'center of mass') 
    stimCoM = All(ind).out.exp.stimCoM;
    numCells = size(All(ind).out.exp.zdfData,1);
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
% else
%     StimDistance = All(ind).out.anal.StimDistance;
% end

%% Find the distance in feature space between cells and targets
% Presented angles
oriVals = [NaN 0:45:315];

% All possible targets and their orientation preferences
allTargets = All(ind).out.exp.targetedCells;
targetedCellOris = nan(size(allTargets));
targetedCellOris(~isnan(allTargets)) = oriVals(All(ind).out.anal.prefOri(allTargets(~isnan(allTargets))));

% All cells in the FOV and their orientation preferences
cellOris = oriVals(All(ind).out.anal.prefOri);

% Feature distance between cells and targets (in radians)
ThetaDistance = abs(targetedCellOris-cellOris);
ThetaDistance(ThetaDistance>180) = abs(ThetaDistance(ThetaDistance>180)-360);
ThetaDistance(ThetaDistance==135)=45;
ThetaDistance(ThetaDistance==180)=0;

% If either a target or cell prefer the gray screen (denoted by nan, 
% assign a distance of 90 degrees
ThetaDistance(isnan(ThetaDistance))=90;

ThetaDistance = ThetaDistance*pi/180;

%%

%load in some variables you'll need
numStims = numel(All(ind).out.exp.stimParams.Seq);
offTargetRisk = All(ind).out.anal.offTargetRisk;
ROIinArtifact = All(ind).out.anal.ROIinArtifact;
respMat = All(ind).out.anal.respMat;
baseMat = All(ind).out.anal.baseMat;


%     now iterate through every stim and see the response as a function of
%     distance
i= ensID;
holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble

if i~=1
    
    Tg=All(ind).out.exp.rois{holo};
    dists = StimDistance(Tg,:);
    thetaDists = ThetaDistance(Tg,:);
        
    %             minDistbyHolo(i,:) = minDist;
    %             geoDistbyHolo(i,:) = geoDist;
    %             meanDistbyHolo(i,:) = meanDist;
    %             harmDistbyHolo(i,:) = harmDist;
    
    switch distType
        case 'min'
            minDist = min(dists,[],1);
            distToUse = minDist;
        case 'geo'
            try
                geoDist = geo_mean(dists,1);
            catch
                geoDist = geomean(dists,1);
            end
            distToUse = geoDist;
        case 'mean'
            meanDist = mean(dists,1);
            distToUse = meanDist;
        case 'harm'
            harmDist = harmmean(dists,1);
            distToUse = harmDist;
        case 'center of mass'
            center_of_mass = mean(stimLoc(All(ind).out.exp.rois{holo},:));
            centerDist = sqrt(sum((center_of_mass-allLoc).^2,2))';
            distToUse = centerDist;
        case 'conn dist'
            
            sigma_dist_rossi = 142; % in microns
            r_0 = 0.44;
            r_p = 0.6730;
            sigma_theta_rossi = 0.3665; % in radians
            
            cr = (r_0 + r_p *exp(-thetaDists.^2/(2*sigma_theta_rossi^2)));
            expDist = exp(-dists.^2/(2*sigma_dist_rossi^2));
            conndist = sum(cr.*expDist);
            distToUse = conndist;
            
        case 'conn dist v2'
            
            sigma_dist_rossi = 142; % in microns
            sigma_dist_narrow = 25;
            kappa_narrow = 0.075;
            r_0 = 0.44;
            r_p = 0.6730;
            sigma_theta_rossi = 0.3665; % in radians
            
            cr = (r_0 + r_p *exp(-thetaDists.^2/(2*sigma_theta_rossi^2)));
%             cr = 1;
            expDist = kappa_narrow*exp(-dists.^2/(2*sigma_dist_narrow^2))+(1-kappa_narrow)*exp(-dists.^2/(2*sigma_dist_rossi^2));
            conndist = sum(cr.*expDist);
            distToUse = conndist;
            
        case 'gauss'
            sigma_dist_rossi = 142;
            conndist = sum(exp(-dists.^2/(2*sigma_dist_rossi^2)),1);
            distToUse = conndist;
        otherwise
            disp('dist type not understood, using min')
            minDist = min(dists,[],1);
            distToUse =minDist;
    end
    
    for d = 1:numel(distBins)-1
        cellsToUse = ~ROIinArtifact' &...
            ~offTargetRisk(holo,:) &...
            distToUse > distBins(d) &...
            distToUse <= distBins(d+1) &...
            cellToUseLimit;
        %if you eventually want to use this with different vis
        %stims use the commented out code.
        %                 popRespDistEns(c,v,d) = nanmean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
    
        popRespDistEns(d) = nanmean(squeeze(respMat(i,1,cellsToUse) - baseMat(i,1,cellsToUse)));
    end
    

    cellsToUse = ~ROIinArtifact' &...
        ~offTargetRisk(holo,:) &...
        cellToUseLimit;
    %if you eventually want to use this with different vis
    %stims use the commented out code.
    %                 popRespDistEns(c,v,d) = nanmean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
    
    distVals = distToUse(cellsToUse);
    
    cellResp = squeeze(respMat(i,1,cellsToUse) - baseMat(i,1,cellsToUse));
    
else
    popRespDistEns = nan([numel(distBins)-1 1]);
end


%     popRespDistEns{ind}=popRespDist;

if plotFlag
    figure(4);clf
    p =plot(popRespDistEns);
    p.LineWidth=2;
    ylim([-0.3 +0.3])
end
