function [dat, All] = pullDat(All,depth)
if ~isfield(All.out,['dat' num2str(depth)])...
        || isempty(eval(['All.out.dat' num2str(depth)])) 
    
    disp('Loading Dat...') 
    
    try
    temp = load(fullfile(All.out.info.path,...
        ['F_' All.out.info.mouse '_' All.out.info.date '_plane' num2str(depth) '_proc.mat']));
    dat = temp.dat;
    eval(['All.out.dat' num2str(depth) ' = dat;']);
    catch
        disp('Error loading')
        disp(fullfile(All.out.info.path,...
        ['F_' All.out.info.mouse '_' All.out.info.date '_plane' num2str(depth) '_proc.mat']))
        dat=[];
    end
    
    
else
%     disp('dat Found Returning')
    dat = eval(['All.out.dat' num2str(depth)]);
end