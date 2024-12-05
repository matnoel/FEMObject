function s = calc_depla(S,q,varargin)
% function s = calc_depla(S,q,varargin)
% S : MODEL
% q : deplacement

if israndom(S) || israndom(q)
    error(' calc_deplapc pas programmé')
else
    q = unfreevector(S,q);
    
    [s,fieldstorage,fieldtype] = calc_elemfield(S,@depla,q,varargin{:});
    
end