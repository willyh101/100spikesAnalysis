%%
% compare max vs circular mean for determining PO
po = outVars.prefOris{1};
oris = 0:45:315;

osi = outVars.osi{1};
osi(po==1)=[];
po(po==1)=[];

podeg = oris(po-1);

figure(2)
clf
subplot(1,2,1)
histogram(podeg, 9)
title('PO by Max')

subplot(1,2,2)
histogram(outVars.circTuning{1}, 9)
title('PO by Circular Mean')
%%

figure(1)
clf
subplot(1,3,1)
scatter(podeg, outVars.circTuning{1}, [], osi, 'filled')
ylabel('Circular PO')
xlabel('Max PO')
c = colorbar;
c.Label.String = 'OSI';

subplot(1,3,2)
scatter(podeg, outVars.circTuning{1}, [], outVars.circVar{1}, 'filled')
ylabel('Circular PO')
xlabel('Max PO')
c = colorbar;
c.Label.String = 'Circ. Var (deg)';

% do OSI on new fancy tuning method
prefOri = outVars.circTuning{1};
ortho1 = mod(prefOri - 90, 135);
ortho2 = mod(prefOri + 90, 135);
orthoOri = cat(1,ortho1, ortho2);
% but have to cast back to 0:45:135
prefOriBinned = interp1(oris, oris, prefOri, 'nearest', 'extrap');
orthoOriBinned = interp1(oris, oris, orthoOri, 'nearest', 'extrap');
% then go back to visID to get idxs
for i=1:numel(prefOriBinned)
    o = prefOriBinned(i);
    oo = orthoOriBinned(:,i);
    prefIDs(i) = find(o==oris);
    orthoIDs1(i) = find(oo(1,:)==oris);
    orthoIDs2(i) = find(oo(2,:)==oris);
end
orthoIDs = cat(1,orthoIDs1, orthoIDs2);

%%
oriCurveBL = outVars.circCurves{1}';
% oriCurveBL = curves - min(curves);
OSI=[];
for i=1:numel(prefOriBinned)
    OSI(i) = (oriCurveBL(prefIDs(i),i) - mean(oriCurveBL(orthoIDs(:,i)',i)))...
        / (oriCurveBL(prefIDs(i),i) + mean(oriCurveBL(orthoIDs(:,i)',i)));
end

subplot(1,3,3)
scatter(podeg, outVars.circTuning{1}, [], OSI, 'filled')
ylabel('Circular PO')
xlabel('Max PO')
c = colorbar;
c.Label.String = 'OSI by circ tuning';