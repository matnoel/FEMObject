function P = setdefaultparam(P,varargin)
% function P = setdefaultparam(P,'paramname1',paramvalue1,'paramname2',paramvalue2,...)
% rajoute les parametres s'il ne sont pas deja presents

paramnames = varargin(1:2:end-1);
paramvalues = varargin(2:2:end);

for i=1:length(paramnames)
    if ~isfield(P.param,paramnames{i})
        P.param = setfield(P.param,paramnames{i},paramvalues{i});
    end
end
