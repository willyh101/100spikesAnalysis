function [cost] = pred_dist_effect(D,R,reg,X);
    X = X./reg; %regularizing change 
    
        Eeffect = sum(gaussian(D,X(1),X(2)));
        Ieffect = sum(gaussian(D,X(3),X(4)));
        
        NetEffect = Eeffect-Ieffect;
        cost = double(sum((R' - NetEffect).^2));
end