function [p, Rsq] = simplifiedLinearRegression(A,B)
%Ian was annoyed at how many lines of code were neded to say that something
%was linearly related so functionalized it.


p  = polyfit(A,B,1);
f = polyval(p,A);
Bbar = mean(B);
SStot = sum((B - Bbar).^2);
SSreg = sum((f - Bbar).^2);
SSres = sum((B - f).^2);
R2 = 1 - SSres/SStot;
R = corrcoef(A,B);
Rsq = R(1,2).^2;