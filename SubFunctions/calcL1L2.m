function [L1 L2 L3] =  calcL1L2(testData,ExpectedData,method)
%return the L1 and L2 values for a dataset size Cells x Obs
if nargin <3
    method =1;
end
meanResp = mean(ExpectedData,2);

if method==1
Err = meanResp - testData ;
else
Err = testData; %dan says this should be included, so you look at the
end
% population sparsity instead of the sparsity of the change
Err = mean(Err'); 

L1 = sum(abs(Err));
L2 = sqrt(sum(Err.^2));
L3 = sum(abs(Err).^3).^(1/3);
