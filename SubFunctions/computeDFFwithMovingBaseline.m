function [dfData, zdfData] =  computeDFFwithMovingBaseline(allData,smFactor,interpRange,percAsBaseline)
% calculate df/f and zscored df/f from Ca Data. 
% allData is a cells x frames or cells x frames x trials of calcium
% fluorescence
% smFactor is the span of the moving baseline

% Written by Ian Oldenburg 2019. Adesnik Lab. UC Berkeley

if nargin<2
    smFactor =1000; %baseline size
    interpRange = 10; % to save time, only calulate baseline every n frames and interp
    percAsBaseline = 10; % percentile to use ase baseline value
elseif nargin<3
    interpRange = 10; % to save time, only calulate baseline every n frames and interp
    percAsBaseline = 10; % percentile to use ase baseline value
elseif nargin<4
    percAsBaseline = 10; % percentile to use ase baseline value
end

%unroll allData
sz = size(allData);
if numel(sz)==3
allDataUnroll = reshape(allData,[sz(1) sz(2)*sz(3)]);
totFrame = sz(2)*sz(3);

elseif numel(sz)==2
    allDataUnroll = allData;
    totFrame = sz(2);
else
    disp('Error Wrong dimmensionality')
end
    
numCells = sz(1);


% calculate the moving baseline 
t = tic;
framesToDo = [1:interpRange:totFrame-1 totFrame];
mvPData = nan([numCells numel(framesToDo)]);
i=0;
parfor k=1:numel(framesToDo);
    i=framesToDo(k);
    % determine window to baseline
    win = round([max((i-smFactor/2),1) : min(i+smFactor/2,totFrame)]);
    % take baseline as nth percentile of each cell in that window
    mvPData(:,k) = prctile(allDataUnroll(:,win),percAsBaseline,2);
    
end

% interperlate baseline to full size
f0Data = zeros(size(allDataUnroll));
for i=1:numCells
    f0Data(i,:) = interp1(framesToDo,mvPData(i,:),1:totFrame,'linear');
end

%correct negative numbers and compute dF/F
minF0Data = min(f0Data(:)) - 1;

allDataToUse = allDataUnroll - minF0Data;
f0DataToUse = f0Data - minF0Data;
dfDataUnroll = (allDataToUse-f0DataToUse)./f0DataToUse;

%Z-Score
zdfDataUnroll = zscore(dfDataUnroll,[],2);

% reshape to original size and output
dfData = reshape(dfDataUnroll,sz);
zdfData = reshape(zdfDataUnroll,sz);



disp(['done. Took ' num2str(toc(t)) ' s'])