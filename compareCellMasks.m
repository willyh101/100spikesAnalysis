% figure(1414);clf
% [ro co]= splitPretty(numExps,4,3)

for ind = 4; 1:numExps
%     subplot(ro,co,ind)
    
    
    
 
    cellCat1 = find(All(ind).out.info.defNeg);
    cellCat2 = find(All(ind).out.info.defPos);
    
    figure(2);clf  
        ax = subplot(1,1,1);
    ImDat = zeros([15 15 numel(cellCat1)]);
    
    imSize = 25;
    ImDatRot1 = zeros([imSize imSize numel(cellCat1)]);
    for i=1:numel(cellCat1);
      
                fprintf('.')

        [All(ind), outIm] = plotCellMask(All(ind),cellCat1(i),ax,1);
        sz = size(outIm);
        ImDat(1:sz(1),1:sz(2),i)=outIm;
        
        stats = regionprops(outIm>0,{'Orientation' 'Eccentricity'});
        EccsCell1(i) =  stats.Eccentricity;
        OriCell1(i) = stats.Orientation;

        
        Or = stats.Orientation;
        
        newIM = imrotate(outIm,Or);
        sz = size(newIM);
        newIM = padarray(newIM,[floor((imSize-sz(1))/2) floor((imSize-sz(2))/2)]);
        sz = size(newIM);
        ImDatRot1(1:sz(1),1:sz(2),i) = newIM;
        

    end
    figure(6);imagesc(mean(ImDat,3))
    title('Opsin Negative')
        figure(8);imagesc(mean(ImDatRot1,3))
    title('Opsin Negative')
    
    disp('Now Opsin Positive')
    
     ImDat2 = zeros([15 15 numel(cellCat2)]);
         ImDatRot2 = zeros([imSize imSize numel(cellCat1)]);

    for i=1:numel(cellCat2);
        fprintf('.')

        [All(ind), outIm] = plotCellMask(All(ind),cellCat2(i),ax,1);
        sz = size(outIm);
        ImDat2(1:sz(1),1:sz(2),i)=outIm;
        
        stats = regionprops(outIm>0,{'Orientation' 'Eccentricity'});
        EccsCell2(i) =  stats.Eccentricity;
        OriCell2(i) = stats.Orientation;
        
                newIM = imrotate(outIm,Or);
        sz = size(newIM);
        newIM = padarray(newIM,[floor((imSize-sz(1))/2) floor((imSize-sz(2))/2)]);
        sz = size(newIM);
        ImDatRot2(1:sz(1),1:sz(2),i) = newIM;
        
    end
    figure(7);imagesc(mean(ImDat2,3))
    title('Opsin Positive')
            figure(9);imagesc(mean(ImDatRot2,3))
    title('Opsin Positive')
    
    disp('done')
    
end