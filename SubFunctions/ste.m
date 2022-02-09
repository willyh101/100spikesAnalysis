function e = ste(v)

e = nanstd(v(:))./sqrt(numel(v));
