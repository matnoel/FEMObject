function levels = LEVELSETS(varargin)
% function levels = LEVELSETS(varargin)

if nargin==1 && isa(varargin{1},'LEVELSETS') && ~isa(varargin{1},'LSCRACK')
    levels = varargin{1};
else    
    levels.LS=cell(1,0);
    for k=1:nargin
        if isa(varargin{k},'LEVELSET') || isa(varargin{k},'LSCRACK') 
            levels.LS = [levels.LS , {varargin{k}}];
        elseif isa(varargin{k},'LEVELSETS')
            levels.LS = [levels.LS , varargin{k}.LS];    
        elseif isa(varargin{k},'cell')
            levels_k = LEVELSETS(varargin{k}{:}) ;
            levels.LS = [levels.LS , levels_k.LS]; 
        end
    end
    levels.n = length(levels.LS);
    levels = class(levels,'LEVELSETS');
    levels = setnumber(levels);
    levels = unique(levels);
end
