function ls = LEVELSET(value,varargin)
% function ls = LEVELSET(value,varargin)

if nargin == 1 && isa(value,'LEVELSET')
    ls = value ;
elseif isa(value,'function_handle')
    ls.value = [{value} , {varargin}] ;
    ls.D = [] ;
    ls.number=[];
    ls.sign = 1;
    ls.material = [];
    ls.enrich = 0 ;
    ls.enrichtype = 0;
    ls=class(ls,'LEVELSET');
else
    ls.value = value ;
    ls.D = getclassin('MODEL',varargin);
    if ~isa(value,'cell') && isempty(ls.D)
        % error('rentrer un MODEL pour definir une levelset')
    end
    
    if isa(ls.value,'MULTIMATRIX') && length(ls.value)==1 && ~israndom(ls.value)
        ls.value = double(ls.value);
    end
    
    ls.number = [];
    ls.sign = 1;
    ls.material = getclassin('MATERIAL',varargin);
    ls.enrich = 0;
    ls.enrichtype = 0;
    
    ls=class(ls,'LEVELSET');
end
