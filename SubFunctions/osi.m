function result = osi(tuning_curve)

[~, prefOri]= max(tuning_curve);

orthoOri = prefOri-2;
orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;

orthoOri2 = orthoOri+4;
orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;

orthoOri = cat(1,orthoOri, orthoOri2);

oriCurveBL = oriCurve - min(oriCurve);

OSI=[];
for i=1:numel(prefOri)
    OSI(i) = (oriCurveBL(prefOri(i),i) - mean(oriCurveBL(orthoOri(:,i)',i)) )...
        / (oriCurveBL(prefOri(i),i)+ mean(oriCurveBL(orthoOri(:,i)',i)) );
    OSI(prefOri==1)=nan;
end

result = OSI;