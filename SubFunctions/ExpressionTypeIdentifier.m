function [All,outVars] = ExpressionTypeIdentifier(All,outVars)
 numExps = numel(All);
 
 mouseNameList = {
     'I127_1'
     'HB40_1'
     'HB47'
     'I132_2'
     'I135_1'
     'I136_1'
     'W4_3'
     'I137_3'
     'I138_1'
     'W20_1'
     'mora_tre'
     'w14_1'
     'w21_1'
     'i139_2'
     'HB95'
     };
 
 ExpressionTypeList = {
     'PHP CAG'
     'AAV Syn'
     'IUE CAG'
     'AAV CamK2'
     'AAV CamK2'
     'AAV Syn'
     'AAV Syn'
     'AAV Syn'
     'AAV Syn'
     'Ai203'
     'AAV Tre'
     'AAV Tre'
     'PHP Tre' 
     'AAV Tre'
     'Ai203'
     };
     
 [uniqueExpressionTypes, ~, ExpressionTypeNum] = unique(ExpressionTypeList);
     
indExpressionType = zeros([1 numExps]);
 for ind =1:numExps
     mouseName = All(ind).out.info.mouse;
     
     idx = find(cellfun(@(x) strcmp(lower(mouseName),lower(x)),mouseNameList));
     if isempty(idx);
         disp(['Error mouse: ' mouseName '. Not Detected...']);
     else
         All(ind).out.info.ExpressionType = ExpressionTypeList{idx};
         All(ind).out.info.ExpressionTypeNum = ExpressionTypeNum(idx);
         indExpressionType(ind) =  ExpressionTypeNum(idx);
     end
 end
 
 outVars.uniqueExpressionTypes = uniqueExpressionTypes;
 outVars.indExpressionType=indExpressionType;