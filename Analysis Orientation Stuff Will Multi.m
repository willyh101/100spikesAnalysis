%% Load Experiments

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');

numExps = numel(loadList);

clear All
for ind = 1:numExps
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
end

%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=6;
recWinSec = [1.25 2.5];

 for ind =1:numExps
     All(ind).out.anal.numCells = size(All(ind).out.exp.zdfData,1);
     numCells(ind) = size(All(ind).out.exp.zdfData,1);
     
     if ~isfield(All(ind).out.info,'FR')
         All(ind).out.info.FR=FRDefault;
     end
     
     sz = size(All(ind).out.exp.zdfData);
     
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
     rdata = squeeze(mean(All(ind).out.exp.zdfData(:,winToUse,:),2));
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);
     bdata = squeeze(mean(All(ind).out.exp.zdfData(:,bwinToUse,:),2));
     
     All(ind).out.exp.rdData=rdata;
     All(ind).out.exp.bdata=bdata;

     sz2 = size(All(ind).out.vis.zdfData);
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz2(2)]) ;
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);

     rdata = squeeze(mean(All(ind).out.vis.zdfData(:,winToUse,:),2));
     bdata = squeeze(mean(All(ind).out.vis.zdfData(:,bwinToUse,:),2));
     
     All(ind).out.vis.rdata=rdata;
     All(ind).out.vis.bdata=bdata;
     
     temp = unique([All(ind).out.exp.holoTargets{:}]);
     temp(isnan(temp))=[];
     All(ind).out.anal.targets = temp;
     numUniqueTargets(ind) =numel(temp);
 end
%% Determine the OSI from the Vis section of each cell.

for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    uVisID = unique(All(ind).out.vis.visID);
    uVisID(uVisID==0)=[];
    
    oriCurve=[];
    for i=1:numel(uVisID)
        v= uVisID(i);
        trialsToUse = All(ind).out.vis.visID==v & All(ind).out.vis.lowMotionTrials;
        
        oriCurve(i,:)=mean(All(ind).out.vis.rdata(:,trialsToUse),2);
    end
    
    All(ind).out.anal.oriCurve = oriCurve;
    
    [maxOriVal maxOriIndex]= max(oriCurve);
    All(ind).out.anal.prefOri = maxOriIndex;
    
    % (Rpref ? Rorth)/(Rpref + Rorth)
    % visID 1 is always catch (I Hope... Will is more confident...ish)
    
    prefOri = maxOriIndex;
    orthoOri = prefOri-2;
    orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
    
    orthoOri2 = orthoOri+4;
    orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
    
    orthoOri = cat(1,orthoOri, orthoOri2);
    
    
    
    %     orthoOri(prefOri==1)=NaN;
    
    oriCurveBL = oriCurve - min(oriCurve);%oriCurve(1,:);
    
    OSI=[];
    for i=1:numel(prefOri)
        OSI(i) = (oriCurveBL(prefOri(i),i) - mean(oriCurveBL(orthoOri(:,i)',i)) ) / (oriCurveBL(prefOri(i),i)+ mean(oriCurveBL(orthoOri(:,i)',i)) );
        %     OSI = (oriCurveBL(prefOri) - oriCurveBL(orthoOri) ) ./ ( oriCurveBL(prefOri)+oriCurveBL(orthoOri) )
        OSI(prefOri==1)=nan;
    end
    
    All(ind).out.anal.OSI=OSI;
    
    
    pVisR=[];pVisT=[];
    for i=1:All(ind).out.anal.numCells
        trialsToUse = All(ind).out.vis.visID~=0 & All(ind).out.vis.lowMotionTrials;
        pVisR(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
        
        trialsToUse = All(ind).out.vis.visID~=0 & All(ind).out.vis.visID~=1 & All(ind).out.vis.lowMotionTrials;
        pVisT(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
    end
    
    All(ind).out.anal.pVisR = pVisR;
    All(ind).out.anal.pVisT = pVisT;
    
    alpha = 0.05;
    
    All(ind).out.anal.visPercent = sum(pVisR<alpha) / numel(pVisR);
    visPercent(ind) =  All(ind).out.anal.visPercent;

    meanOSI=[];ensembleOSI=[];ensembleOriCurve =[];ensemblePref=[];
    
    for i=1:numel(All(ind).out.exp.holoTargets)
        ht = All(ind).out.exp.holoTargets{i};
        ht(isnan(ht))=[];
        meanOSI(i)=nanmean(OSI(ht));
        
        ensOriCurve = mean(oriCurve(:,ht),2);
        ensembleOriCurve(i,:)= ensOriCurve;
        [maxOriVal maxOriIndex]= max(ensOriCurve);
        prefOri = maxOriIndex;
        
        orthoOri = prefOri-2;
        orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
        orthoOri2 = orthoOri+4;
        orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
        orthoOri = cat(1,orthoOri, orthoOri2);
        oriCurveBL = ensOriCurve - min(ensOriCurve);
        
        ensembleOSI(i) = (oriCurveBL(prefOri)- mean(oriCurveBL(orthoOri'))) / (oriCurveBL(prefOri)+ mean(oriCurveBL(orthoOri')));
        ensemblePref(i) = prefOri;
    end
    
    deg = [nan 0:45:315];
    All(ind).out.anal.ensembleOSI=ensembleOSI;
    All(ind).out.anal.meanOSI=meanOSI;
    All(ind).out.anal.ensembleOSI=ensembleOSI;
    All(ind).out.anal.ensembleOriCurve=ensembleOriCurve;
    All(ind).out.anal.ensemblePrefOri=ensemblePref;
    All(ind).out.anal.ensemblePrefDeg=deg(ensemblePref);
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% Pretty plots of OSI and tunings

clear allOSI ensOSI meanOSI
% OSI across all cells, all experiments
for i = 1:numel(All)
    allOSI{i} = All(i).out.anal.OSI(:);
    ensOSI{i} = All(i).out.anal.ensembleOSI(:);
    meanOSI{i} = All(i).out.anal.meanOSI(:);
end

% unroll
allOSI = cell2mat(allOSI(:));
ensOSI = cell2mat(ensOSI(:));
meanOSI = cell2mat(meanOSI(:));

% plot
figure(1)
clf(1)
colors = {rgb('royalblue'), rgb('firebrick'), rgb('coral')};

hold on
h(1) = histogram(allOSI, 25);
h(2) = histogram(ensOSI, 25);
h(3) = histogram(meanOSI, 25);

for i=1:numel(h)
    h(i).Normalization = 'pdf';
    h(i).BinWidth = 0.05;
    h(i).FaceColor = colors{i};
    h(i).FaceAlpha = 0.44;
    kde(i) = fitdist(h(i).Data, 'kernel');
    p = plot(h(i).BinEdges, pdf(kde(i),h(i).BinEdges));
    p.LineWidth=2;
    p.Color= colors{i};
end

hold off

set(gcf(),'Name','OSI Across All Expts')
title('OSI Across All Expts')
xlabel('Orientation Selectivity Index')
ylabel('PDF')
legend('All Cells', 'Ensemble', 'Ensemble Mean')

%% now look at tuned ensembles
% get plots for the 2 different methods with higher bin count
f2 = figure(2);
clf(f2)
hold on
h(1) = histogram(ensOSI, 50);
%h(2) = histogram(meanOSI, 50);

% try using ensembleOSI > 0.3 for "tuned", this is arbitrary
OSIthreshold = 0.3;
for i = 1:numel(All)
    istuned = All(i).out.anal.ensembleOSI > OSIthreshold;
    tunedEnsembles = All(i).out.exp.holoTargets(istuned);
    untunedEnsembles = All(i).out.exp.holoTargets(~istuned);
    tunedEnsembleIdx = find(All(i).out.anal.ensembleOSI => OSIthreshold);
    untunedEnsembleIdx = find(All(i).out.anal.ensembleOSI < OSIthreshold);
    All(i).out.anal.tunedEnsembles = tunedEnsembles;
    All(i).out.anal.untunedEnsembles = untunedEnsembles;
    All(i).out.anal.tunedEnsembleIdx = tunedEnsembleIdx;
    All(i).out.anal.untunedEnsembleIdx = untunedEnsembleIdx;
end

% plot the OSIs of tuned vs untuned ensembles
for i = 1:numel(All)
    allOSI{i} = All(i).out.anal.OSI(:);
    ensOSI{i} = All(i).out.anal.ensembleOSI(:);
    meanOSI{i} = All(i).out.anal.meanOSI(:);
end

% unroll
allOSI = cell2mat(allOSI(:));
ensOSI = cell2mat(ensOSI(:));
meanOSI = cell2mat(meanOSI(:));

f3 = subplots(1,2,1)

% and thier preferred oris

    
    














