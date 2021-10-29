

ensToUse = find(outVars.ensemblesToUse & outVars.ensOSI>0.7 & outVars.meanEnsOSI>0.5);
% 
  ensToUse = find(outVars.ensemblesToUse & outVars.ensOSI<0.3 & outVars.meanEnsOSI>0.5);
 ensToUse = find(outVars.ensemblesToUse & outVars.ensOSI<0.3 & outVars.meanEnsOSI<0.5);
ensToUse = [536 518 614];

figure(19);clf

for i =1:numel(ensToUse)
    subplot(1,3,i)
    
    hold on
    e = ensToUse(i);
    ind = ensIndNumber(e);
    
    s = outVars.ensHNumber(e);
    h = All(ind).out.exp.stimParams.roi{s};
    %     All(ind).out.exp.stimParams
    htg = All(ind).out.exp.holoTargets{h};
    htg(isnan(htg))=[];
    tuninCurves = outVars.tuningCurves{ind}(:,htg);
    tuninCurveSEM = outVars.tuningCurvesSEM{ind}(:,htg);
    
    clist = fliplr(colorMapPicker(5,'viridis'));
    ors=[];
    for k = 1:size(tuninCurves,2)
        %        eb =  errorbar(2:9,tuninCurves(2:9,k),tuninCurveSEM(2:9,k));
        eb =  plot(1:9,tuninCurves(1:9,k));
        
        prefOri =outVars.prefOris{ind}(htg(k));
        prefOri(prefOri>5) = prefOri-4;
        ors(k)=prefOri;
        %        eb.CapSize=0;
        eb.LineWidth = 1;
        %        eb.Color = [rgb('grey') 0.5];
        eb.Color = [clist{prefOri} 0.5];
    end
    ors
    
    eb = errorbar(2:9,outVars.ensCurve(2:9,e),outVars.ensCurveSEM(2:9,e));
    eb.CapSize =0;
    eb.LineWidth = 2;
    eb.Color = 'r';
    
    eb = errorbar(1,outVars.ensCurve(1,e),outVars.ensCurveSEM(1,e));
    eb.CapSize =0;
    eb.LineWidth = 2;
    eb.Color = 'k';
    
    xticks([1 3 5 7 9])
    xticklabels({'c' '0' '90' '180' '270'})
    

    title(e)
    ylim([-0.5 1.5])
%    pause 
end
outVars.ensOSI(ensToUse)
outVars.meanEnsOSI(ensToUse)
