ensemblesToUseList = find(outVars.ensemblesToUse...
    & outVars.ensOSI > 0.75...
    & outVars.meanEnsOSI > 0.5...
    );

disp([num2str(numel(ensemblesToUseList)) ' Ensembles Included'])

for E = 1:numel(ensemblesToUseList);
    try
        ens = ensemblesToUseList(E);
        ind = outVars.ensIndNumber(ens);
        hNum = outVars.ensHNumber(ens);
        ht = All(ind).out.exp.holoTargets{hNum};
        ht(isnan(ht)) = [];
        %
        figure(237);clf

        numTargets = outVars.numCellsEachEns(ens);
        [wi hi] = splitPretty(numTargets,3,4);
        for k = 1:numel(ht)
            h = ht(k);

            subplot(hi,wi,k)

            %     errorbar(outVars.tuningCurves{ind}(:,h),outVars.tuningCurvesSEM{ind}(:,h));
            errorbar(All(ind).out.anal.oriCurve(:,h),All(ind).out.anal.oriCurveSEM(:,h));
            cellOSI = outVars.osi{ind}(h);
            title(['Cell : ' num2str(h) '. osi ' num2str(cellOSI,2)]);
        end

        figure(238);clf

        hi = 1; wi = 4;
        subplot(hi,wi,1)
        e = errorbar(outVars.ensCurve(:,ens),outVars.ensCurveSEM(:,ens));
        xticks([1:9]);
        xticklabels([-10 0:45:315])
        xtickangle(45)

        e.LineWidth = 2;
        ensOSI = outVars.ensOSI(ens);
        title(['Composite Tuning Curve, Ensemble : ' num2str(ens) '. osi ' num2str(ensOSI,2)]);
        xlabel('Ori')
        ylabel('Average Zdf')


        subplot(hi,wi,2)

        allCellsMeanResp = outVars.mRespEns{ens};
        exclCells = All(ind).out.anal.offTargetRisk(hNum,:) | All(ind).out.anal.ROIinArtifact';
        allCellsMeanResp(exclCells) = nan;
        plot(sort(allCellsMeanResp),'.')
        refline(0)
        title('Mean Response of all non stim cells')
        xlabel('cell ID')
        ylabel('\Delta zdf')
        xlim([0 sum(~isnan(allCellsMeanResp))])

        subplot(hi,wi,3)
        respOSI = outVars.osi{ind};
        [sRespOSI sidx] = sort(respOSI);
        plot(sRespOSI, allCellsMeanResp(sidx),'.')
        xlabel('Responder OSI')
        ylabel('\Delta zdf')
        refline(0)
        title(['Mean Resp: ' num2str(nanmean(allCellsMeanResp),2) '. signrank: ' num2str(signrank(allCellsMeanResp),2)]);

        subplot(hi,wi,4)
        respOri = outVars.prefOris{ind};
        [sRespOri sidxori] = sort(respOri);
        plot(sRespOri, allCellsMeanResp(sidxori),'.')
        xlabel('Responder prefOri')
        ylabel('\Delta zdf')
        refline(0)
        xticks([1:9]);
        xticklabels([-10 0:45:315])
        xtickangle(45)

        pause
        
    catch
        disp(['Something broken with ens number ' num2str(ens) ' (experiment ' num2str(ind) ', stim ' num2str(hNum) '.'])
        disp('Proceeding to next ensemble...')
    end
end

