function out = makePaddedTS(data, startFrame, alignTo, finalLength)

newStart = startFrame-alignTo+1;

dat = squeeze(mean(data,3));
dat = mean(dat,1);

if startFrame < alignTo
    padby = alignTo - startFrame;
    pdat = padarray(dat', padby, 'pre');
    pdat(pdat == 0) = nan;
    dat = pdat';
    newStart = newStart + padby;
end

if finalLength > length(dat)
    padby = finalLength - length(dat);
    pdat = padarray(dat', padby, 'post');
    pdat(pdat == 0) = nan;
    dat = pdat';
end

dat = dat(newStart:finalLength);

out = dat - mean(dat(1:5),2);