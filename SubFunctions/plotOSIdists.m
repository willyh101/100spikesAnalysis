function plotOSIdists(outVars, opts)

method = opts.ensOSImethod;

figure(43)
clf
colors = {rgb('royalblue'), rgb('firebrick'), rgb('coral')};


allOSI = cell2mat(outVars.osi(:)');
allOSI = allOSI(cell2mat(outVars.isVisR(:)'));

hold on

hg(1) = histogram(allOSI, 25);

try
    hg(2) = histogram(outVars.(method), 25);
catch
    error('Not a valid ensemble osi calculation method.')
end


for i=1:numel(hg)
    hg(i).Normalization = 'pdf';
    hg(i).BinWidth = 0.05;
    hg(i).FaceColor = colors{i};
    hg(i).FaceAlpha = 0.44;
    kde(i) = fitdist(hg(i).Data', 'kernel');
    p = plot(hg(i).BinEdges, pdf(kde(i),hg(i).BinEdges));
    p.LineWidth=2;
    p.Color= colors{i};
end

hold off

set(gcf(),'Name','OSI Across All Expts')
title('OSI Across All Expts')
xlabel('Orientation Selectivity Index')
ylabel('PDF')
legend('All Cells', 'Ensemble OSI')