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
     'i140_2'
     'w26_1'
     'w18_3'
     'w19_1'
     'w29_1'
     'w29_3'
     'w30_2'
     'w32_2'
     'I143'
     'I147'
     'i151_3'
     'w37_2'
     'I154_2'
     'I154_1'
     'I156_1'
     'HB42_1'
     'W40_2'
     'I149_2'
     'I158_1'
     'I162_1'
     'I162_2'
     'W42_2'
     'W46_2'
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
     'AAV Tre'
     'Ai203'
     'AAV Tre'
     'AAV Tre'
     'AAV Tre'
     'AAV Tre 2s'
     'AAV Tre 2s'
     'AAV Tre 2s'
     'Ai203'
     'Ai203'
     'AAV Tre 2s'
     'neo-IV Tre 2s'
     'AAV Tre 2s'
     'AAV Tre 2s'
     'SepW1 CAG 2s'
     'control'
     'AAV Tre 2s'
     'control'
     'SepW1 CAG 2s'
     'SepW1 CAG 2s'
     'SepW1 CAG 2s'
     'AAV Tre 2s'
     'AAV Tre 2s'
     };

  genoTypeNameList = {
     'double'
     'EMX'
     'iue'
     'PV'
     'PV'
     'EMX'
     'EMX'
     'EMX'
     'EMX'
     'Ai203'
     'double'
     'SST'
     'double'
     'PV'
     'Ai203'
     'PV'
     'Ai203'
     'SST'
     'SST'
     'SST'
     'SST'
     'PV'
     'SST'
     'Ai203'
     'Ai203'
     'PV'
     'EMX'
     'double'
     'double'
     'SepW1'
     'Control'
     'SST'
     'Control'
     'SepW1'
     'SepW1'
     'SepW1'
     'SST'
     'SST'
     };
 
 [uniqueExpressionTypes, ~, ExpressionTypeNum] = unique(ExpressionTypeList);

  [uniqueGenoTypes, ~, GenoTypeNum] = unique(genoTypeNameList);

 
indExpressionType = zeros([1 numExps]);
indGenoType = zeros([1 numExps]);
 for ind =1:numExps
     mouseName = All(ind).out.info.mouse;

     idx = find(cellfun(@(x) strcmp(lower(mouseName),lower(x)),mouseNameList));
     if isempty(idx);
         disp(['Error mouse: ' mouseName '. Not Detected...']);
     else
         All(ind).out.info.ExpressionType = ExpressionTypeList{idx};
         All(ind).out.info.ExpressionTypeNum = ExpressionTypeNum(idx);
         indExpressionType(ind) =  ExpressionTypeNum(idx);
         
          All(ind).out.info.GenoType = genoTypeNameList{idx};
         All(ind).out.info.GenoTypeNum = GenoTypeNum(idx);
         indGenoType(ind) =  GenoTypeNum(idx);
         
     end
 end

 outVars.uniqueExpressionTypes = uniqueExpressionTypes;
 outVars.indExpressionType=indExpressionType;

  outVars.uniqueGenoTypes = uniqueGenoTypes;
 outVars.indGenoType=indGenoType;
 