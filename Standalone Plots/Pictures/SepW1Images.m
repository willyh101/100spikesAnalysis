[fn pn] = uigetfile;
%%

dat = ScanImageTiffReader(fullfile(pn,fn)).data;


G = dat(:,:,1:2:end);
R = dat(:,:,2:2:end); 

nDepths =3;


%%
depth = 2;

Gi = G(:,:,depth:nDepths:end);
Ri = R(:,:,depth:nDepths:end);

GiBackup = Gi;
RiBackup = Ri;

%%
Gi = GiBackup;
Ri = RiBackup;


img = zeros([512 512 3]);
gRange = [50 225 2];
rRange = [30 50 5];

Gi= double(Gi);
Gi = Gi-gRange(1); 
Gi(Gi<0)=0;
Gi = Gi./(gRange(2)-gRange(1));
Gi(Gi>1)=1;
Gi = Gi.*gRange(3);

Ri = double(Ri);
Ri = Ri-rRange(1);
Ri(Ri<0)=0;

Ri = Ri./(rRange(2)-rRange(1));
Ri(Ri>1)=1;
Ri = Ri.*rRange(3);

img(:,:,2) = mean(Gi,3);
img(:,:,1) = mean(Ri,3);

figure(1);clf

subplot(2,3,[1 2 4 5])
image(img)
axis square
box off
axis off

hold on
width = 100 /800*512;
r = rectangle('position',[20 490 width 10]); 
r.FaceColor = rgb('lightgrey');
r.LineStyle = 'none';

subRange = 15;
loc1 = [80 304];

subplot(2,3,3)
subIm = img(round(loc1(1)-subRange/2):round(loc1(1)+subRange/2),...
    round(loc1(2)-subRange/2):round(loc1(2)+subRange/2), :);
image(subIm)
axis square
axis off

width = 10 /800*512;
r = rectangle('position',[2 14 width 1]); 
r.FaceColor = rgb('lightgrey');
r.LineStyle = 'none';


loc1 = [149 149];%[172 261 ];
subplot(2,3,6)
subIm = img(round(loc1(1)-subRange/2):round(loc1(1)+subRange/2),...
    round(loc1(2)-subRange/2):round(loc1(2)+subRange/2), :);
image(subIm)
axis square
axis off
width = 10 /800*512;
r = rectangle('position',[2 14 width 1]); 
r.FaceColor = rgb('lightgrey');
r.LineStyle = 'none';
