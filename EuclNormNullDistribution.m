tempResp = rand([10000 8]);


for i=1:size(tempResp,1)
        vector = tempResp(i,1:8);
        vector = vector-min(vector);
        vector = vector./max(vector);
        EN = sqrt(sum((vector.^2)));
        
        EN = 1 - ((EN-1)./(sqrt(8)-1));
        
        EucNormAll(i) = EN;
end
    
figure(26);clf
histogram(EucNormAll,100)
xlim([0 1])