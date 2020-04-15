function s = sem2(data, dim)

if nargin < 2
    dim = 1;
end

s = nanstd(data, [], dim)./(sqrt(sum(~isnan(data), dim)));