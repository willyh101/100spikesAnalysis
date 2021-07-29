IdsToUse = find(outVars.ensemblesToUse);

for II = 1:numel(IdsToUse)
ensID = IdsToUse(II);

ind = outVars.ensIndNumber(ensID);

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

%load in some variables you'll need
numStims = numel(All(ind).out.exp.stimParams.Seq);
offTargetRisk = All(ind).out.anal.offTargetRisk;
ROIinArtifact = All(ind).out.anal.ROIinArtifact;
respMat = All(ind).out.anal.respMat;
baseMat = All(ind).out.anal.baseMat;


%     now iterate through every stim and see the response as a function of
%     distance
H = outVars.ensHNumber(ensID);
% for i= 1:numStims
holo = All(ind).out.exp.stimParams.roi{H}; % Better Identifying ensemble

if H~=1
    
    Tg=All(ind).out.exp.rois{holo};
    dists = StimDistance(Tg,:);
else
    dists=[];
end

if isempty(CellToUseVar)
    cellToUseLimit = ones([1 size(All(ind).out.anal.respMat,3)]);
elseif islogical(CellToUseVar)
    cellToUseLimit=CellToUseVar;
else
    try
        cellToUseLimit = eval(['All(ind).out.' CellToUseVar]);
    catch
        disp('ERROR variable not found')
        %         break
    end
end

cellsToUse = ~ROIinArtifact' &...
    ~offTargetRisk(holo,:) &...
    cellToUseLimit;


D = dists(:,cellsToUse);
R = squeeze(respMat(H,1,cellsToUse) - baseMat(H,1,cellsToUse));

reg = [10 1 10 1]; %regularization change
fun = @(x) pred_dist_effect(D,R,reg,x);
x0 = [1 50 1 50];
x0 = x0.*reg;
%x order EAmp ESignma IAmp ISigma

opts.Display = 'off';% 'iter';
opts.DiffMinChange = 1e-3;
% opts.CheckGradients = 'false';
% [x,fval,exitflag,output,grad,hessian] = fminunc(fun,x0,opts); 


A = [];
b = [];
Aeq = [];
beq = [];

lb = [0,0,0,0];
ub = [inf,inf,inf,inf];

x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],opts);
x=x./reg; 
disp(['Ens ' num2str(II) ' of ' num2str(numel(IdsToUse)) ' : ' num2str(x,3)]);
GT(:,II) = x;
end

%%
figure(10);clf
subplot(1,2,1)
plot(GT(2,:))
hold on
plot(GT(4,:));
title('Sigma')

subplot(1,2,2)
plot(GT(1,:))
hold on
plot(GT(3,:));
title('Amplitude')

median(GT(2,:))
median(GT(4,:))
signrank(GT(2,:),GT(4,:))


median(GT(1,:))
median(GT(2,:))
signrank(GT(1,:),GT(2,:))



