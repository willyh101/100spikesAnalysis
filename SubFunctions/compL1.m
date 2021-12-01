function L1 = compL1(A)
A(isnan(A))=[];
L1 = sum(abs(A))/numel(A);