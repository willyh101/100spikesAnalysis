function kdehist(data, binwidth, color)

    h = histogram(data);
    h.Normalization = 'pdf';
    h.BinWidth = binwidth;
    h.FaceAlpha = 0.4;
    h.FaceColor = color;
    kde = fitdist(h.Data, 'kernel');
    p = plot(h.BinEdges, pdf(kde,h.BinEdges));
    p.LineWidth = 2;
    p.Color = color;
