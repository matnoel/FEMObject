function [indicator_U,indicator_w,indicator_lambda] = initializeIndicators(m,patches)
% function [indicator_U,indicator_w,indicator_lambda] = initializeIndicators(m,patches)
% Initializes indicators
% m: maximum number of parameters (iterations or samples or ...)
% patches: Patches

n = numel(patches);

indicator_U = zeros(1,m);
indicator_w = repmat({zeros(1,m)},1,n);
indicator_lambda = repmat({zeros(1,m)},1,n);

end
