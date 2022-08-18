ensToUse = find(outVars.ensemblesToUse);



% for i  = 1:numel(ensToUse)
%     e = ensToUse(i);
%     ind = ensIndNumber(e);
%     
%     s = outVars.ensHNumber(e); %i think this is the order of resp
%     h = All(ind).out.exp.stimParams.roi{s};
%     htg = All(ind).out.exp.holoTargets{h};
%     htg(isnan(htg))=[];
%     
%     rtg = All(ind).out.exp.rois{h};
%     
%     us = unique(All(ind).out.exp.stimID);
% 
%     sID = us(s); %StimID value for this ensemble
%     
%     squeeze(All(ind).out.anal.resp(s,:,:))
%     
%     
%     
% end
%%
plotIt=0;

clear allNegRate allPosRate
allNegRate =[];
allPosRate =[];
for ii = 1:numel(IndsUsed)
    ind = IndsUsed(ii);
    
    theseEns = find(ensemblesToUse & ensIndNumber==ind);
    
    clear holoResp holoSTE
    for i  = 1:numel(theseEns)
        e = theseEns(i);
        ind = ensIndNumber(e);
        
        s = outVars.ensHNumber(e); %i think this is the order of resp
        h = All(ind).out.exp.stimParams.roi{s};
        htg = All(ind).out.exp.holoTargets{h};
        htg(isnan(htg))=[];
        
        rtg = All(ind).out.exp.rois{h};
        
        us = unique(All(ind).out.exp.stimID);
        
        sID = us(s); %StimID value for this ensemble
        
        temp = squeeze(All(ind).out.anal.respMat(s,:,:));
        temp(outVars.offTargetRiskEns{e} | All(ind).out.anal.ROIinArtifact')=nan;
        temp2 = squeeze(All(ind).out.anal.baseMat(s,:,:));
        
        temp = temp - temp2;
        
        holoResp(i,:) =temp;
        
        temp = squeeze(All(ind).out.anal.stdMat(s,:,:));
        temp(outVars.offTargetRiskEns{e} | All(ind).out.anal.ROIinArtifact')=nan;
        
        temp2 =  squeeze(All(ind).out.anal.numMat(s,:,:));
        temp = temp./sqrt(temp2); 
        holoSTE(i,:) =temp;
        
        isNeg2 = double(holoResp<0 & holoSTE*2.58+holoResp<0);
        isPos2 = double(holoResp>0 & holoResp-holoSTE*2.58>0);
        
        isNeg2(isnan(holoResp))=NaN; 
        isPos2(isnan(holoResp))=NaN;

        
        isNonSigNeg = holoResp<0 & holoSTE*1.96+holoResp>0;
        isNonSigPos = holoResp>0 & holoResp-holoSTE*1.96<0;
    end

    allNegRate = cat(1,allNegRate,nanmean(isNeg2,2));
        allPosRate = cat(1,allPosRate,nanmean(isPos2,2));

    
        if plotIt
        
    for i=18; %1:size(holoResp,2)
        cellID = i;
        
        figure(22);clf
        e = errorbar(holoResp(:,cellID),holoSTE(:,cellID)*2.58);%*1.96);
        title(cellID);
        xlabel('Hologram Number')
        ylabel('Evoked dF/F')
        xlim([0.5 size(holoResp,1)+0.5])        
        refline(0)
        
        e.LineWidth=1;
        e.Color=rgb('grey');
        e.Marker='o';
        e.LineStyle='none';
        
        pause
    end
    
    for i= 1:size(holoResp,1)
        figure(23);clf
        H = i;
        
        thisHoloResp = holoResp(H,:);
        thisHoloCI = holoSTE(H,:)*2.58;%*1.96;
        
        [sortHoloResp sidx] = sort(thisHoloResp);
        sortHoloCI = thisHoloCI(sidx);
        
        e = errorbar(sortHoloResp,sortHoloCI);
        e.LineWidth=1;
        e.Color=rgb('grey');
        e.Marker='none';
        e.LineStyle='none';
        e.CapSize=0;
        
        isNeg = sortHoloResp<0 & sortHoloCI+sortHoloResp<0;
        isPos = sortHoloResp>0 & sortHoloResp-sortHoloCI>0;
        
        hold on
        e = errorbar(find(isNeg | isPos),sortHoloResp(isNeg | isPos),sortHoloCI(isNeg | isPos));
        e.Color=rgb('red');
        e.Marker='none';
        e.LineStyle='none';
        e.CapSize=0;
        
        sum(isNeg)
        sum(isPos)
        pause
    end
        end

end


figure(4);clf
hold on
h1 = histogram(allNegRate,40);
h2 = histogram(allPosRate,40);
