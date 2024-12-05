function ls = enrich(ls,type,varargin)
% function ls = enrich(ls,type,varargin)

if any(type>0)
ls.enrich = 1;
ls.enrichtype = type;
end
