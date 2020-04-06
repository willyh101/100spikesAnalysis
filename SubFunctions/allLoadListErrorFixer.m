function [All] = allLoadListErrorFixer(All)
%Hard Code Error Fixes

%note the prefered way of catching errors
nameToUse = '190826_I132_2_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    All(indToUse).out.exp.visIDBackup=All(indToUse).out.exp.visID;
    visID = All(indToUse).out.exp.visIDBackup;
    visID(visID==2)=6;
    All(indToUse).out.exp.visID=visID;
end

%% has all 0 (grey screen)
nameToUse = '200217_w14_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    All(indToUse).out.exp.visID = ones(1,length(All(indToUse).out.exp.visID));
end

%% experiment xx somehow has a 30 and 3 in the begining of stimParams
nameToUse = '200224_w14_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    badID = All(indToUse).out.exp.stimID==16;
    All(indToUse).out.exp.stimID(badID)=[];
    All(indToUse).out.exp.zdfData(:,:,badID)=[];
    All(indToUse).out.exp.allData(:,:,badID)=[];
    All(indToUse).out.exp.runVal(badID,:)=[];
    All(indToUse).out.exp.lowMotionTrials(badID)=[];
    All(indToUse).out.exp.visID(badID)=[];
    All(indToUse).out.exp.outputsInfo.OutputStims(2)=[];
    All(indToUse).out.exp.outputsInfo.OutputNames(2)=[];
    All(indToUse).out.exp.outputsInfo.OutputOrder(2)=[];
    All(indToUse).out.exp.outputsInfo.OutputPatterns(2)=[];
    All(indToUse).out.exp.stimParams.Seq(2)=[];
    All(indToUse).out.exp.stimParams.numPulse(2)=[];
    All(indToUse).out.exp.stimParams.roi(2)=[];
end

%% Copied from Ori 
%% fix Exp1 two unused holos
nameToUse = '190418_I127_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    All(indToUse).out.exp.stimParams.Seq([3 4])=[];
    All(indToUse).out.exp.stimParams.numPulse([3 4])=[];
    All(indToUse).out.exp.stimParams.roi([3 4])=[];
    All(indToUse).out.exp.stimParams.Hz([2 3])=[];
    All(indToUse).out.exp.stimParams.numCells([2 3])=[];
    All(indToUse).out.exp.holoTargets([2 3])=[];
end
%% fix Exp10; ignoring intermediate contrasts
nameToUse = '191017_I136_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
     All(indToUse).out.exp.visIDold = All(indToUse).out.exp.visID;
    visID = All(indToUse).out.exp.visID;
    newVisID = visID;
    newVisID(visID==0)=1;
    newVisID(visID==1)=2;
    newVisID(visID~=0 & visID~=1)=0;
    All(indToUse).out.exp.visID = newVisID;
end

%% fix Exp15; merge two blank visIDs
nameToUse = '191217_mora_tre_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    All(indToUse).out.exp.visIDold = All(indToUse).out.exp.visID;
    visID = All(indToUse).out.exp.visID;
    newVisID = visID;
    newVisID(visID==4)=1;
    All(indToUse).out.exp.visID = newVisID;
    
    All(indToUse).out.info.offsets = [1.4667 -88.8667];

end

%% Eliminate low contrast Exp 16
nameToUse = '191230_w20_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    All(indToUse).out.exp.visID(All(indToUse).out.exp.visID(1,:)~=1)=0;
    All(indToUse).out.info.offsets = [-0.2 -0.1333];
end

%% fix Exp 9 visID 
nameToUse = '190606_HB47_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
All(indToUse).out.exp.visID = [ones([1 425]) ones([1 400])*2];
All(indToUse).out.exp.visCond = cat(2,repmat ([1;nan],[1 425]),All(indToUse).out.exp.visCond);
end


