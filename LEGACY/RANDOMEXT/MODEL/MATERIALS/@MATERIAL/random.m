function mat = random(mat,varargin)
% function mat = random(mat,varargin)

RV = RANDVARS(mat);
RVsample = random(RV,varargin{:});

mat = randomeval(mat,RVsample,RV);
