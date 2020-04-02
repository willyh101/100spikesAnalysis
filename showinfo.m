function t = showinfo(All)

for i=1:numel(All)
    mouse{i} = All(i).out.info.mouse;
    date{i} = All(i).out.info.date;
    t1{i} = All(i).out.info.epochText1;
    t2{i} = All(i).out.info.epochText2;
    n{i} = i;
end

c={n{:}; mouse{:}; date{:}; t1{:}; t2{:}}';
t = cell2table(c, 'VariableNames', {'Number', 'Mouse', 'Date', 'Epoch Text 1', 'Epoch Text 2'});