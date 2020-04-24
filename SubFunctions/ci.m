function CI = ci(data, dim)

if nargin < 2
    dim = 1;
end

sem = std(data, dim)/sqrt(size(data, dim));
m_data = mean(data, dim);
ts = tinv([0.025, 0.975], size(data, dim) - 1);
CI = 