function curve = makeTuningCurve180(tuning_curve)
% send the tuning curve WITHOUT the blank condition
% must come in format of: condition x cell

assert(size(tuning_curve,1)==8, ['Number of conditions must be 8, not ' ...
    num2str(size(tuning_curve,1)) '. Remove any catch conditions and ensure '...
    'tuning curve is condition/orientation by cell.'])

data_temp(:,:,1) = tuning_curve(1:4,:);
data_temp(:,:,2) = tuning_curve(5:8,:);
curve = mean(data_temp,3);