function [All, outVars] = getEnsembleOSI(All, outVars)

clear meanEnsOSI meanEnsOri cmeanEnsOri ensCurve ensCurveSEM ensOSI ensPO

for ind = 1:numel(All)
    
    clear meanOSI meanOri cmeanOri curve2save curveSEM curve ensPO ensOSI

    for h = 1:numel(All(ind).out.exp.holoTargets)
        % ht = holo targets
        ht = All(ind).out.exp.holoTargets{h};
        ht(isnan(ht)) = [];
        
        % get mean OSI, easy
        meanOSI(h) = nanmean(All(ind).out.anal.osi(ht));
        
        % get mean Ori by simple mean
        meanOri(h) = nanmean(idx2ori(All(ind).out.anal.prefOri(ht), [nan 0:45:315]));
        
        % mean ori by circular mean
        pos = idx2ori(All(ind).out.anal.prefOri(ht), [nan 0:45:315]);
        pos(isnan(pos)) = [];
        radPref = circ_ang2rad(pos);
        cmeanOri(h) = circ_rad2ang(circ_axial(circ_mean(radPref')));
        
        % mean ori by ensemble tuning curve
        curve = mean(All(ind).out.anal.oriCurve(:, ht), 2);
        curve2save(:,h) = curve;
        curveSEM(:, h) = sem2(All(ind).out.anal.oriCurve(:, ht), 2);
        [~, ensPOidx] = max(curve);
        curve(:, ensPOidx==1) = nan;
        curve(1, :) = [];
        ensPO(h) = idx2ori(ensPOidx, [nan 0:45:315]);
        ensOSI(h) = osi(curve);
    end
    
    All(ind).out.exp.meanEnsOSI = meanOSI;
    All(ind).out.exp.meanEnsOri = meanOri;
    All(ind).out.exp.cmeanEnsOri = cmeanOri;
    All(ind).out.exp.ensCurve = curve2save;
    All(ind).out.exp.ensCurveSEM = curveSEM;
    All(ind).out.exp.ensOSI = ensOSI;
    All(ind).out.exp.ensPO = ensPO;
    
    
    meanEnsOSItemp{ind} = meanOSI;
    meanEnsOritemp{ind} = meanOri;
    cmeanEnsOritemp{ind} = cmeanOri;
    ensCurvetemp{ind} = curve2save;
    ensCurveSEMtemp{ind} = curveSEM;
    ensOSItemp{ind} = ensOSI;
    ensPOtemp{ind} = ensPO;
end

outVars.meanEnsOSI = cell2mat(meanEnsOSItemp(:)');
outVars.meanEnsOri = cell2mat(meanEnsOritemp(:)');
outVars.cmeanEnsOri = cell2mat(cmeanEnsOritemp(:)');
outVars.ensCurve = cell2mat(ensCurvetemp(:)');
outVars.ensCurveSEM = cell2mat(ensCurveSEMtemp(:)');
outVars.ensOSI = cell2mat(ensOSItemp(:)');
outVars.ensPO = cell2mat(ensPOtemp(:)');