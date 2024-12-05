function ls = actualise(ls,varargin)
% function ls = actualise(ls,varargin)

if ~iseval(ls) && isclassin('MODEL',varargin)
    ls = lseval(ls,getclassin('MODEL',varargin));
elseif iseval(ls) && isclassin('MODEL',varargin)
    % warning('on ne stocke plus le model dans la levelset')
    % ls.D = getclassin('MODEL',varargin);
elseif ~iseval(ls)
    error('utiliser lseval sur un MODEL')
end

