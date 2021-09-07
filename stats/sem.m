function s = sem(vector)

s = nanstd(vector)/(sqrt(sum(~isnan(vector))));