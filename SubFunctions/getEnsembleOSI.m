function [All, outVars] = getEnsembleOSI(All, outVars)

%Get the number of spikes in each stimulus

for ind = 1:numel(All)

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
    
    
    outVars.meanEnsOSI{ind} = meanOSI;
    outVars.meanEnsOri{ind} = meanOri;
    outVars.cmeanEnsOri{ind} = cmeanOri;
    outVars.ensCurve{ind} = curve2save;
    outVars.ensCurveSEM{ind} = curveSEM;
    outVars.ensOSI{ind} = ensOSI;
    outVars.ensPO{ind} = ensPO;
end