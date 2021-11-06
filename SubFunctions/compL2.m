function L2 = compL2(A)
A(isnan(A))=[];
L2 = sqrt(sum(A.^2))/numel(A);
