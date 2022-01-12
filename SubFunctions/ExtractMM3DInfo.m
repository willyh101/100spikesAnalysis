
useMM3Dim=1;
    
    rootpath = dat.ops.RootDir;
    dataPath = fullfile(rootpath,'Other');
    
    fList = dir(dataPath);
    
    if useMM3Dim
        
        filePointer=[];
        for i =1:numel(fList)
            a = strfind(fList(i).name,'makeMasks3D_img');
            if ~isempty(a) %naturally will find last 1020 image
                filePointer =i;
            end
        end
        
        if isempty(filePointer)
            if ~isdir(dataPath)
                a = strsplit(dataPath,'\');
                dataPath = fullfile(a{1:end-1});
            end
            [fn fp] = uigetfile([dataPath '\*.mat']);
            temp = load(fullfile(fp,fn));
            
        else
            temp = load(fullfile(fList(filePointer).folder,fList(filePointer).name));
        end
        clear opsinImage
        
        for i = 1:nDepthsTotal
            try
            imToUse2 = temp.img{i}(:,:,1)'; %red channel
            catch
                imToUse2 = temp.img1{i}(:,:,1)'; %red channel
            end
            opsinImage{i} = imToUse2;
        end
        out.info.OpsinImage = opsinImage;
    else
        filePointer=[];
        for i =1:numel(fList)
            a = strfind(fList(i).name,'_1020_');
            if ~isempty(a) %naturally will find last 1020 image
                filePointer =i;
            end
        end
        
        disp('Loading Red Stack')
        
        if isempty(filePointer)
            if ~isdir(dataPath)
                a = strsplit(dataPath,'\');
                dataPath = fullfile(a{1:end-1});
            end
            [fn fp] = uigetfile([dataPath '\*.tif']);
            raw1020Img = ScanImageTiffReader(fullfile(fp,fn)).data;
        else
            raw1020Img = ScanImageTiffReader(fullfile(fList(filePointer).folder,fList(filePointer).name)).data;
            
        end
        raw1020Img=raw1020Img(:,:,2:2:end);
        
        disp('Aligning')
        nDepthsTotal = 3; %eventually don't want to hardcode it
        for i = 1:nDepthsTotal
            disp(['Depth: ' num2str(i) ' ']);
            imToUse = raw1020Img(:,:,i:nDepthsTotal:end);
            
            %                     mn = prctile(imToUse(:),0.5);%min(imToUse(:));% prctile(imToUse(:),0.5);
            %                     mx = prctile(imToUse(:),100);%max(imToUse(:))-100;%prctile(imToUse(:),99.5);
            %
            %                     imToUse = (imToUse-mn)./(mx-mn);
            
            [imToUse2, dxs, dys] = simpleAlignTimeSeries (imToUse); %align to account for any motion
            %                                     imToUse = mean(imToUse,3);
            imToUse2 = mean(imToUse2,3);
            opsinImage{i} = imToUse2;
        end
    end
    out.info.OpsinImage = opsinImage;
    
    disp('done.')
    %%extract 1020 value
    value1020=[];value1020_2=[];opsinNegative=[];
    [gridX gridY] = ndgrid(1:size(imToUse2,1),1:size(imToUse2,2));
    
    radiusToUse =3;
    
    figure(1);clf
    for i = 1:nDepthsTotal
        subplot(1,nDepthsTotal,i)
        imagesc(opsinImage{i})
        
        try
            offsetsToUse = out.exp.offsets;
        catch
            offsetsToUse = All(ind).out.info.offsets;
        end
        
        
        dCoM = out.exp.allCoM(out.exp.allDepth==i,:)-offsetsToUse;
        thisPlanevals=[];
        for k=1:size(dCoM,1)
            circle(dCoM(k,1),dCoM(k,2),radiusToUse,rgb('lime'));
            
            X = gridX-dCoM(k,2);
            Y = gridY-dCoM(k,1);
            L = sqrt(X.^2+Y.^2)<=radiusToUse;
            
            value1020(end+1)= mean(opsinImage{i}(L));
            thisPlanevals(k) = mean(opsinImage{i}(L));
        end
        
        %         threshold = prctile(opsinImage{i}(:),75);
        threshold = max(thisPlanevals)/100;
        opsinNegative =  [opsinNegative thisPlanevals<threshold];
        
        
        %             unqDepth = unique(out.exp.holoRequest.targets(:,3));
        %             dCoM = out.exp.holoRequest.targets(out.exp.holoRequest.targets(:,3)==unqDepth(i),1:2);
        dCoM = out.exp.stimCoM(out.exp.stimDepth==i,:)-offsetsToUse;
        
        for k=1:size(dCoM,1)
            circle(dCoM(k,1),dCoM(k,2),radiusToUse,rgb('red'));
            
            X = gridX-dCoM(k,2);
            Y = gridY-dCoM(k,1);
            L = sqrt(X.^2+Y.^2)<=radiusToUse;
            
            value1020_2(end+1)= mean(opsinImage{i}(L));
        end
        drawnow
    end
    disp('done')
    %         threshold = 10; % value to call def negative
    
    figure(3);clf
    subplot(1,2,1)
    hold on
    plot(sort(value1020),'.');
    tCell = unique(out.exp.targetedCells);
    tCell(isnan(tCell))=[];
    
    plot(value1020(tCell),'.')
    
    plot(value1020_2,'.')
    
    line([0 numel(value1020)],[threshold threshold]);
    % set(gca,'yscale','log')
    
    subplot(1,2,2)
    h = histogram(value1020,1000);
    
    %     xline(threshold);
    line([threshold threshold],[0  max(h.Values)]);
    
    out.info.value1020=value1020;
    out.info.value1020Targets = value1020_2;
    
    %     out.info.opsinNegative = value1020<threshold;
    out.info.opsinNegative = opsinNegative;
    
    disp([ num2str(sum(out.info.opsinNegative)) ' definitely opsin Negative Cells. ' num2str(mean(out.info.opsinNegative)*100) '%']);
    
    disp('press enter')