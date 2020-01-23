%% identify the vis stim conditions...

% start by displaying all the epoch text so I know what's there
cellswithdata = [1 4 6 7 8 9]
for i=cellswithdata%1:numel(All)
    disp(['Expt: ' num2str(i) ', mouse: ' All(i).out.info.mouse ', date: ' All(i).out.info.date '...'])
    disp(All(i).out.info.epochText1')
    disp(All(i).out.info.epochText2')
    %disp(All(i).out.exp.vis_oris)
    %disp(All(i).out.exp.vis_contrasts)
    pause 
    disp('---')
end

%% manual error checking
for i=1:numel(All)
    disp('----')
    disp(['Expt: ' num2str(i) ', mouse: ' All(i).out.info.mouse ', date: ' All(i).out.info.date '...'])
    uniqueVis = unique(All(i).out.exp.visID);
    disp(['unique vises: ' num2str(uniqueVis)])
    in = input('does this seem right?? [y/n]  ', 's');
    if in == 'y'
        oris_in = input('ok what are the oris?  ', 's');
        All(i).out.exp.vis_oris = oris_in;
        contrasts_in = input('what about contrasts?  ', 's');
        All(i).out.exp.vis_contrasts = contrasts_in;
    end
end
%% put the stim conditions into the correct order/format
% should be contrast L->H, ori L->H

%cellswithdata = [2 3 5 12 14 15 16]
cellswithdata = [3 5 12 14 15 16];

for i=cellswithdata
    All(i).out.exp.vis_oris = str2num( All(i).out.exp.vis_oris);
    All(i).out.exp.vis_contrasts = str2num(All(i).out.exp.vis_contrasts);
end

%%
%clear vis_out
cellswithdata = [9];
for i=cellswithdata
    vis_key = [];
    oris = 90;%sort(All(i).out.exp.vis_oris, 'MissingPlacement','First');
    cons = 1;%sort(All(i).out.exp.vis_contrasts);
    c=1;
    for j=1:numel(cons)
        for k=1:numel(oris)
            vis_key(:,c) = [cons(j); oris(k)];
            c=c+1;
        end
    end
    vis_out{i} = vis_key;
end