function plotEnsPO(outVars)


figure(45)
clf
colors = {rgb('royalblue'), rgb('firebrick'), rgb('coral')};
oris = [nan 0:45:315];

allPO = idx2ori(cell2mat(outVars.prefOris(:)'), oris);
% allPO = allPO(cell2mat(outVars.isVisR(:)'));

hold on


subplot(2,2,1)
histogram(allPO, 9)
title('All Cells Pref Ori')
xlabel('Preferred Orientation')
ylabel('Count')

subplot(2,2,2)
histogram(outVars.meanEnsOri, 9)
title('Ensemble Cells Mean PO')
xlabel('Preferred Orientation')
ylabel('Count')

subplot(2,2,3)
histogram(outVars.cmeanEnsOri, 9)
title('Ensemble Cells Circular Mean PO')
xlabel('Preferred Orientation')
ylabel('Count')

subplot(2,2,4)
histogram(outVars.ensPO, 9)
title('Mean Ensemble PO')
xlabel('Preferred Orientation')
ylabel('Count')