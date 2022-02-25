function out = fwhm(c1)

% divide c1 by sqrt(2) bc matlab reports it multiplied by the by sqrt(2) to
% get actual sigma
out = 2*sqrt(2*log(2))*c1/sqrt(2);