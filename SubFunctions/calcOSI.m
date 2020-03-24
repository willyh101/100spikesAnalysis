function [All, outVars] = calcOSI(All, outVars)

% settings
numExps = numel(All);

if ~isfield(All(1).out.anal, 'oriCurve')
    error("No ori curve, can't calculate OSI.")
end
    

% method 1, the "old" way
for ind = 1:numExps
    % get this expts variables for All struct
    oriCurve = All(ind).out.anal.oriCurve;
    prefOri = All(ind).out.anal.prefOri;
    orthoOri = All(ind).out.anal.orthoOri;
    numCells = All(ind).out.anal.numCells;
    
    % no negative numbers please
    oriCurveBL = oriCurve - min(oriCurve);
    
    % calculate OSI
    osi = [];
    for i = 1:numCells
        osi(i) = (oriCurveBL(prefOri(i),i) - mean(oriCurveBL(orthoOri(:,i)',i))) ...
            / (oriCurveBL(prefOri(i),i) + mean(oriCurveBL(orthoOri(:,i)',i)) );
        osi(prefOri==1) = nan;
    end
    
    All(ind).out.anal.osi = osi;
    osis{ind} = osi;
end

outVars.osi = osis;
