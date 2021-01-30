x1 = outVars.ensOSI';
x1 = x1(ensemblesToUse);
x2 = outVars.ensMaxD';
x2 = x2(ensemblesToUse);
% x3 = outVars.
% y = outVars.popResponseEns(ensemblesToUse);
data = outVars.popResponseDist;
dists = cellfun(@(x) x(:,1), data, 'un', 0);
dists = cell2mat(dists');
y = dists(ensemblesToUse);

X = [ones(size(x1)) x1 x2 x1.*x2];
[b, bint, r, rint, stats] = regress(y,X);

figure(10212)
clf
scatter3(x1,x2,y,'filled')
hold on
x1fit = linspace(min(x1),max(x1),50);
x2fit = linspace(min(x2),max(x2),50);

[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('ensOSI')
ylabel('ensMaxD')
zlabel('popResponseEns')
view(50,10)
hold off

%%


cellDistsEns=[];
cellEnsOSI = [];
cellEnsMaxD = [];
cellEns2Use = [];
cellEnsResp = [];

c = 0;
for ind = 1:numExps
    ROIinArtifact = All(ind).out.anal.ROIinArtifact;
    offTargetRisk = All(ind).out.anal.offTargetRisk;
    pVisR = All(ind).out.anal.pVisR;
    
    
    
    cellDists = All(ind).out.anal.minDistbyHolo;
    
    for j = 1:size(cellDists,1)
        if j ==1
            continue
        end
        c = c+1;
        cellDistsEns{c} = cellDists(j,:);
%         cellDistsEns{c}(offTargetRisk(j-1, :)' | ROIinArtifact) = nan;
        excl = (ROIinArtifact' | offTargetRisk(j-1, :) & pVisR > 0.05);
        
        
        cellEnsOSI{c} = repmat(outVars.ensOSI(c),numel(cellDistsEns{c}),1)';
        cellEnsMaxD{c} = repmat(outVars.ensMaxD(c),numel(cellDistsEns{c}),1)';
        cellEns2Use{c} = repmat(outVars.ensemblesToUse(c), numel(cellDistsEns{c}),1)';
        
        cellDistsEns{c} = cellDistsEns{c}(~excl);
        cellEnsOSI{c} = cellEnsOSI{c}(~excl);
        cellEnsMaxD{c} = cellEnsMaxD{c}(~excl);
        cellEns2Use{c} = cellEns2Use{c}(~excl);
        cellEnsResp{c} = outVars.mRespEns{c}(~excl)';
    	
    end

end

cellDistsEns = cell2mat(cellDistsEns);
cellEnsOSI = cell2mat(cellEnsOSI);
cellEnsMaxD = cell2mat(cellEnsMaxD);
cellEnsResp = cell2mat(cellEnsResp);
cellEns2Use = cell2mat(cellEns2Use);

cellDistsEns = cellDistsEns(cellEns2Use);
cellEnsOSI = cellEnsOSI(cellEns2Use);
cellEnsMaxD = cellEnsMaxD(cellEns2Use);
cellEnsResp = cellEnsResp(cellEns2Use);

%% there are way too many cells so downsample to test shit out

r = randperm(numel(cellEnsResp), 10000);

cellDistsEns = cellDistsEns(r);
cellEnsOSI = cellEnsOSI(r);
cellEnsMaxD = cellEnsMaxD(r);
cellEnsResp = cellEnsResp(r);

dat = table(cellEnsResp', cellDistsEns', cellEnsOSI', cellEnsMaxD');
dat.Properties.VariableNames = {'Response', 'Distance from Target', 'Ensemble OSI', 'Ensemble Spread'};



%%

% x1 = cellDistsEns';
x1 = cellEnsOSI';
x2 = cellEnsMaxD';
y = cellEnsResp';

X = [ones(size(x1)) x1 x2 x1.*x2];
[b, bint, r, rint, stats] = regress(y,X);

figure(10212)
clf
% scatter3(x1,x2,y,'filled')
hold on
x1fit = linspace(min(x1),max(x1),50);
x2fit = linspace(min(x2),max(x2),50);

[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('ensOSI')
ylabel('ensMaxD')
zlabel('popResponseEns')
view(50,10)
hold off

%%