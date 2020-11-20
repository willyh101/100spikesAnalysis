function meanMovie = meanMovieMaker(fileList,targets,color,T,numFrames,depth,nDepths)
%reads and saves all movies specified by targets in fileList
%loading red or green or both and MotionCorrecting 
%T is motion correct from .align file

if nargin<6
    depth =1;
    nDepths =1;
end
% depth = 2;
% nDepths = 2;


%Determine Colors to use
switch color
    case {'red', 'Red', 'R', 'r'}
        G=0;
        R=1;
    case {'green','Green','G','g'}
        G=1;
        R=0;
    case {'Both','both','B','b','RG','GR','rg','gr'}
        G=1;
        R=1;
    otherwise
        G=1;
        R=0;
        fprintf('Did not recognize color command Green only selected...\n');      
end

%cleanup targets to be binary list, checks for all 1 all 0 or list of index
uniqueTargets = unique(targets);
if uniqueTargets == 1;
    fprintf('All Files Selected...\n')
elseif uniqueTargets ==0;
    fprintf('No Files Selected Ending Now.\n');
    meanMovie=[];
    return
elseif numel(uniqueTargets) >2 || ~all(uniqueTargets==[0 1])
    fprintf('Targets not a binary list, converting...')
    newTargets = zeros(max(targets),1);
    for i = 1:max(targets);
        newTargets(i) = any(targets==i);
    end
    targets=newTargets;
end
    
if numel(targets)>numel(fileList)
    fprintf('Too Many Targets, only counting first... \n');
    targets = targets(1:numel(fileList));
elseif numel(targets)<numel(fileList)
    fprintf('Too Few Targets, unspecified files not included... \n');
    padSize = numel(fileList)-numel(targets);
    targets = padarray(targets,[0 padSize],0,'post');
end

fprintf(['Will Extract Movies from ' num2str(sum(targets)) ' files.\n']);
    
%check T
if ~isnumeric(T) || isempty(T)
    fprintf('No MC transform detected. Not MCd\n');
    T=zeros(numel(fileList)*numFrames,2);
end
if length(T)<numel(fileList)*numFrames
    fprintf('MC array not long enough, padded with 0s\n');
    T= padarray(T,numel(fileList)*numFrames-length(T),0,'post');
end


%%
Gout= 0;
Rout=0;

numMOV = sum(targets);

%figure(1);

fprintf(['Extracting ' num2str(numMOV) ' movies:  \n']);
for i=1:numel(fileList)
    if targets(i)
        fprintf([num2str(i) ' ']);
        %[gg, ri]=bigread3(fileList{i});
        v = ScanImageTiffReader(fileList{i}).data();
        v = permute(v,[2 1 3]);
        gg = v(:,:,1:2:end);
        ri = v(:,:,2:2:end);
        
        gg = gg(:,:,depth:nDepths:end);
        if G
            
            ss=size(gg);
            mcg=zeros(ss);
            for n=(1:min(ss(3),numFrames));
                mcg(:,:,n)=circshift(gg(:,:,n),T(n+(numFrames*(i-1)),:));
%                 imagesc(mcg(:,:,n));
%                 pause(0.1);
            end
            if isempty(Gout)
                Gout = double(mcg(:,:,1:numFrames))/numMOV;
            else
                Gout=Gout+double(mcg(:,:,1:numFrames))/numMOV;
            end
        end
        if R
            ss=size(ri);
            mcr=zeros(ss);
            for n=(1:ss(3));
                mcr(:,:,n)=circshift(ri(:,:,n),T(n+(numFrames*(i-1)),:));
            end
            if isempty(Rout)
                Rout = double(mcr)/numMov;
            else
                Rout=Rout+double(mcr)/numMov;
            end
        end
    end
end
 fprintf('\nDone! \n');

if R && G
    meanMovie = cat(5,Gout,Rout);
elseif G
    meanMovie = Gout;
elseif R
    meanMovie = Rout;
else
    meanMovie = Gout;
end

    
