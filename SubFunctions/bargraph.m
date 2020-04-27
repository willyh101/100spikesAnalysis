function b = bargraph(data, varargin)

defaultErrorType = 'sem';
validErrorTypes = {'sem', 'ci', 'std'};
checkErrorType = @(x) any(validatestring(x, validErrorTypes));

p = inputParser;
addOptional(p, 'colors', []); % color specs
addParameter(p, 'errType', defaultErrorType, checkErrorType); % error type selector

parse(p, data, varargin{:})


