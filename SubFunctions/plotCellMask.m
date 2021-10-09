function [All] = plotCellMask(All,cellID,axHandle)

if nargin >2 
    subplot(axHandle)
else
    figure;
end

allDepth =All.out.exp.allDepth;
depth = allDepth(cellID);
cellList = 1:numel(allDepth);
%figure out the ID in as dat sees it
IDinDat = numel(find(allDepth==depth & cellList'<=cellID));

[dat, All] = pullDat(All,depth);

iscelldepth = [dat.stat.iscell];
iscelldepth = find(iscelldepth);
IDinDat = iscelldepth(IDinDat);

sz = size(dat.mimg);
blank= zeros(sz(1:2));

xpix = dat.stat(IDinDat).xpix;
ypix = dat.stat(IDinDat).ypix;
lambda = dat.stat(IDinDat).lam;
npix = dat.stat(IDinDat).npix;

med = dat.stat(IDinDat).med;
xySize = 7;

for i=1:npix
    blank(xpix(i),ypix(i))=lambda(i);
end
sz = size(blank);

blankSmall = blank(max(med(2)-xySize,1):min(med(2)+xySize,sz(1)),...
    max(med(1)-xySize,1):min(med(1)+xySize,sz(2)) );
imagesc(blankSmall)
axis square
axis off
colormap viridis

muPerPx = 800/512;
r = rectangle('Position',[1,xySize*2,10/muPerPx,1]);
r.FaceColor= rgb('white');
r.LineStyle = 'none';

