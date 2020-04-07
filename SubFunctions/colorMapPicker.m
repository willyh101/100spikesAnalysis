function [colorList,thisColor] = colorMapPicker(totalElements,cmap,J);
if nargin<2 || isempty(cmap)
    cmap = colormap('parula');
end
if nargin < 3
    J=1;
end

if isstring(cmap) || ischar(cmap)
    cmap = colormap(cmap);
end

cl=[];
colorList = [];
cl = cmap;
ncs = round(linspace(1,size(cl,1),totalElements));
for i=1:totalElements;
    colorList{end+1} =  cl(ncs(i),:);
end

thisColor=colorList{J};
