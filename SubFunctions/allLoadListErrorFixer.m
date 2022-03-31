function [All] = allLoadListErrorFixer(All,loadList)
%Hard Code Error Fixes

%note the prefered way of catching errors
nameToUse = '190826_I132_2_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    All(indToUse).out.exp.visIDBackup=All(indToUse).out.exp.visID;
    visID = All(indToUse).out.exp.visIDBackup;
    visID(visID==2)=6;
    All(indToUse).out.exp.visID=visID;
end

%% has all 0 (grey screen)
nameToUse = '200217_w14_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    All(indToUse).out.exp.visID = ones(1,length(All(indToUse).out.exp.visID));
end

%% experiment xx somehow has a 30 and 3 in the begining of stimParams
nameToUse = '200224_w14_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    badID = All(indToUse).out.exp.stimID==16;
    All(indToUse).out.exp.stimID(badID)=[];
    All(indToUse).out.exp.zdfData(:,:,badID)=[];
    All(indToUse).out.exp.dfData(:,:,badID)=[];
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
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    All(indToUse).out.exp.stimParams.Seq([3 4])=[];
    All(indToUse).out.exp.stimParams.numPulse([3 4])=[];
    All(indToUse).out.exp.stimParams.roi([3 4])=[];
    All(indToUse).out.exp.stimParams.Hz([2 3])=[];
    All(indToUse).out.exp.stimParams.numCells([2 3])=[];
    All(indToUse).out.exp.rois([2 3])=[];
    All(indToUse).out.exp.holoTargets([2 3])=[];
end

%%
nameToUse = '190415_I127_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    All(indToUse).out.exp.stimParams.Hz=[30 30 30 30];
end

%%
nameToUse = '200224_w14_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    All(indToUse).out.exp.stimParams.Hz(1)=[];
end

%% fix Exp10; ignoring intermediate contrasts
nameToUse = '191017_I136_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
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
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
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
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    All(indToUse).out.exp.visID(All(indToUse).out.exp.visID<4 ...
        & All(indToUse).out.exp.visID>1)=0;
    All(indToUse).out.info.offsets = [-0.2 -0.1333];
end

%% fix Exp 9 visID
nameToUse = '190606_HB47_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    All(indToUse).out.exp.visID = [ones([1 425]) ones([1 400])*2];
    All(indToUse).out.exp.visCond = cat(2,repmat ([1;nan],[1 425]),All(indToUse).out.exp.visCond);
end

%% fix vis sections 190420
nameToUse = '190420_I127_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    sz = size(All(indToUse).out.vis.zdfData);
    All(indToUse).out.vis.runVal = zeros([sz(3) sz(2)]);
    All(indToUse).out.vis.visID = ones([1 sz(3)]);
    disp('Some Vis Data wasn not availble. i guessed')
end

%% fix trial counts 191212 I 138
nameToUse = '191212_I138_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));
if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    sz = size(All(indToUse).out.vis.zdfData);

    All(indToUse).out.vis.runVal = All(indToUse).out.vis.runVal(1:sz(3),:);
    All(indToUse).out.vis.visID = All(indToUse).out.vis.visID(1:sz(3));
end

%% %correct ranging for exp
nameToUse = '200727_W26_1_outfile.mat';
%Not sure if this is right or if it should be 2:65...
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    range = [2:65];
    All(indToUse).out.exp.stimID = All(indToUse).out.exp.stimID(range);
    All(indToUse).out.exp.visID = All(indToUse).out.exp.visID(range);
    All(indToUse).out.exp.runVal = All(indToUse).out.exp.runVal(range,:);
    
    All(indToUse).out.mani.stimID = All(indToUse).out.mani.stimID(range);
    All(indToUse).out.mani.visID = All(indToUse).out.mani.visID(range);
    All(indToUse).out.mani.runVal = All(indToUse).out.mani.runVal(range,:);

    range = [3:44];
    All(indToUse).out.stm.stimID = All(indToUse).out.stm.stimID(range);
        All(indToUse).out.stm.runVal = All(indToUse).out.stm.runVal(range,:);

end



%% merge high and low contrasts in vis epoch
nameToUse = '200728_i140_2_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    All(indToUse).out.vis.visID(All(indToUse).out.vis.visID > 9) ...
        = All(indToUse).out.vis.visID(All(indToUse).out.vis.visID > 9) - 9;
end
%% extra lowRunTrial and visID for some reason
nameToUse = '200807_w26_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    All(indToUse).out.vis.runVal(188,:) = [];
    All(indToUse).out.vis.visID(188) = [];
end
%% fix 'vis' stim to be same as blank holo stim
nameToUse = '210426_w32_2_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    
    if any(All(indToUse).out.exp.stimID == 13)
        All(indToUse).out.exp.stimID(All(indToUse).out.exp.stimID == 13) = 14;
        All(indToUse).out.exp.outputsInfo.OutputStims(1) = [];
        All(indToUse).out.exp.outputsInfo.OutputNames{1} = [];
        All(indToUse).out.exp.outputsInfo.OutputOrder(1) = [];
        All(indToUse).out.exp.outputsInfo.OutputPatterns{1} = [];
        All(indToUse).out.exp.stimParams.Seq(1)=[];
        All(indToUse).out.exp.stimParams.numPulse(1)=[];
        All(indToUse).out.exp.stimParams.roi(1)=[];
        All(indToUse).out.exp.uniqueStims(1) = [];
    end
end

%% fix exp in sang min style awake and anest
nameToUse = '210902_I151_3_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    All(indToUse).out.exp = All(indToUse).out.exp1;
end 

%% fix number of trials error in I143 210426

nameToUse = '210426_I143_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    try All(indToUse).out.vis.runVal = All(indToUse).out.vis.runVal(1:111,:);catch;end
    All(indToUse).out.vis.visID = All(indToUse).out.vis.visID(1:111);
end 

%% fix number of trials error in I154_2 210927

nameToUse = '210927_I154_2_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    disp(['Correcting from Ind: ' num2str(indToUse)]);
    lengthOfDat = size(All(indToUse).out.exp.allData,3);
    try All(indToUse).out.exp.runVal = All(indToUse).out.exp.runVal(2:end,:);catch;end
    All(indToUse).out.exp.visID = All(indToUse).out.exp.visID(2:end);
    All(indToUse).out.exp.stimID = All(indToUse).out.exp.stimID(2:end);
end 


%% fix numCells in 210902_I151_3

nameToUse = '210902_I151_3_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
    All(indToUse).out.exp.stimParams.numCells = All(indToUse).out.exp.stimParams.numCells*10;
end


%% fix vis Window 211014_I156_1

nameToUse = '211014_I156_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
   All(indToUse).out.vis.visStart   = 0.98;
   All(indToUse).out.vis.visStop    = 1.98;
end

%% fix 190521_HB42_1

nameToUse = '190521_HB42_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
  All(indToUse).out.exp.stimParams.numCells = ones([1 60])*10;
end

%% fix 211102

nameToUse = '211102_I158_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));
if ~isempty(indToUse)
    All(indToUse).out.exp2.runVal = All(indToUse).out.exp2.runVal(1:899,:);
    All(indToUse).out.exp2.stimID = All(indToUse).out.exp2.stimID(1:899);
    All(indToUse).out.exp2.visID = All(indToUse).out.exp2.visID(1:899);


end


%% fix vis Window 211014_I156_1

nameToUse = '210517_I147_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
  All(indToUse).out.stm.holoRequest =  All(indToUse).out.stm.holoRequest2;
end

%% 
nameToUse = '210105_W29_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
  All(indToUse).out.stm.runVal =  All(indToUse).out.stm.runVal(1:150,:);
    All(indToUse).out.stm.stimID =  All(indToUse).out.stm.stimID(1:150);

     All(indToUse).out.spk.runVal =  All(indToUse).out.spk.runVal(1:117,:);
    All(indToUse).out.spk.stimID =  All(indToUse).out.spk.stimID(1:117);
        All(indToUse).out.spk.visID =  All(indToUse).out.spk.visID(1:117);

    
end
%%
nameToUse = '200724_W18_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
  All(indToUse).out.exp.stimID = All(indToUse).out.exp.stimID(1:76);
    All(indToUse).out.exp.visID = All(indToUse).out.exp.visID(1:76);
        All(indToUse).out.exp.runVal = All(indToUse).out.exp.runVal(1:76,:);
    All(indToUse).out.exp.lowRunTrials = All(indToUse).out.exp.lowRunTrials(1:76);


end

%%
nameToUse = '210105_W29_1_outfile.mat';
indToUse = find(cellfun(@(x) strcmp(x,nameToUse),loadList));

if ~isempty(indToUse)
  All(indToUse).out.exp.stimID = All(indToUse).out.exp.stimID(1:160);
    All(indToUse).out.exp.visID = All(indToUse).out.exp.visID(1:160);
        All(indToUse).out.exp.runVal = All(indToUse).out.exp.runVal(1:160,:);
    All(indToUse).out.exp.lowRunTrials = All(indToUse).out.exp.lowRunTrials(1:160);


end

%%
disp('Done fixing.')