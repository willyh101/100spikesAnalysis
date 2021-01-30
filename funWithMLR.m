x1 = outVars.ensOSI';
x1 = x1(ensemblesToUse);
x2 = outVars.ensMaxD';
x2 = x2(ensemblesToUse);
% x3 = outVars.
% y = outVars.popResponseEns(ensemblesToUse);
data = outVars.popResponseDist;
dists = cellfun(@(x) x(:,1), data, 'un', 0);
dists = cell2mat(dists');
y = dists(ensemblesToUse);

X = [ones(size(x1)) x1 x2 x1.*x2];
[b, bint, r, rint, stats] = regress(y,X);

figure(10212)
clf
scatter3(x1,x2,y,'filled')
hold on
x1fit = linspace(min(x1),max(x1),50);
x2fit = linspace(min(x2),max(x2),50);

[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('ensOSI')
ylabel('ensMaxD')
zlabel('popResponseEns')
view(50,10)
hold off