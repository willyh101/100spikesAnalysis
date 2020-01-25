function [L1 L2 L3] =  calcL1L2(testData,ExpectedData)
%return the L1 and L2 values for a dataset size Cells x Obs

meanResp = mean(ExpectedData,2);

Err = meanResp - testData ;
% Err = testData;
Err = mean(Err'); 

L1 = sum(abs(Err));
L2 = sqrt(sum(Err.^2));
L3 = sum(abs(Err).^3).^(1/3);
