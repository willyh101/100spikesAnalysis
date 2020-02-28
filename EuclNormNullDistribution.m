tempResp = rand([10140 9]);


for i=1:size(tempResp,1)
        vector = tempResp(i,:);
        vector = vector-min(vector);
        vector = vector./max(vector);
        vector = vector(1:8);
        EN = sqrt(sum((vector.^2)));
        
        EN = 1 - (abs(EN-1)./(sqrt(8)-1));
        
        EucNormNull(i) = EN;
end
    
figure(26);clf
histogram(EucNormNull,100)
% xlim([0 1])