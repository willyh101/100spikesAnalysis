
for i=1:10
    Locs(i,:) = randi(600,[1 2])+[100 100];
end



gridX = 1:2:800;
gridY = 1:2:800;

gridXY=[];
for i = 1:numel(gridY)
    for k = 1:numel(gridX)
        gridXY(:,k,i) = [gridX(k) gridY(i)];
    end
end
sz =size(gridXY);
gridXY = reshape(gridXY,[2 sz(2)*sz(3)]);

distance=[];
for i = 1:10
    fprintf([num2str(i) ' ']);
    for k = 1:size(gridXY,2)
        distance(i,k) = sqrt(sum((Locs(i,:)-gridXY(:,k)').^2));
    end
end

LocsCentroid = mean(Locs);
for k=1:size(gridXY,2)
distToCentroid(k) = sqrt(sum((LocsCentroid-gridXY(:,k)').^2));
end
%%
figure(137);clf

subplot(2,3,1)

minDist = min(distance,[],1); 
scatter(gridXY(1,:),gridXY(2,:),1,minDist')
hold on
s = scatter(Locs(:,1),Locs(:,2));
xlim([1 800])
ylim([1 800])

s.MarkerFaceColor = rgb('red');
colorbar
axis square
title('min')


subplot(2,3,2)
geoDist = geomean(distance,1); 
scatter(gridXY(1,:),gridXY(2,:),1,geoDist')
hold on
s = scatter(Locs(:,1),Locs(:,2));
xlim([1 800])
ylim([1 800])

s.MarkerFaceColor = rgb('red');
colorbar
axis square
title('geometric mean')

subplot(2,3,3)
meanDist = mean(distance,1); 
scatter(gridXY(1,:),gridXY(2,:),1,meanDist')
hold on
s = scatter(Locs(:,1),Locs(:,2));
xlim([1 800])
ylim([1 800])

s.MarkerFaceColor = rgb('red');
colorbar
axis square
title('arithmetic mean')

subplot(2,3,4)
harmDist = harmmean(distance,1); 
scatter(gridXY(1,:),gridXY(2,:),1,harmDist')
hold on
s = scatter(Locs(:,1),Locs(:,2));
xlim([1 800])
ylim([1 800])

s.MarkerFaceColor = rgb('red');
colorbar
axis square
title('harmonic mean')

subplot(2,3,5)
medianDist = median(distance,1); 
scatter(gridXY(1,:),gridXY(2,:),1,medianDist')
hold on
s = scatter(Locs(:,1),Locs(:,2));
xlim([1 800])
ylim([1 800])

s.MarkerFaceColor = rgb('red');
colorbar
axis square
title('median')

subplot(2,3,6)
% EuclDist = sqrt(sum(distance.^2,1)); 
% scatter(gridXY(1,:),gridXY(2,:),1,EuclDist')
% hold on
% s = scatter(Locs(:,1),Locs(:,2));
% xlim([1 800])
% ylim([1 800])
% 
% s.MarkerFaceColor = rgb('red');
% colorbar
% axis square
% title('Euclidean norm')

scatter(gridXY(1,:),gridXY(2,:),1,distToCentroid')
hold on
s = scatter(Locs(:,1),Locs(:,2));
xlim([1 800])
ylim([1 800])

s.MarkerFaceColor = rgb('red');
colorbar
axis square
title('Distance to Centroid')
