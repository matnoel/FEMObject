function mats = MATERIALS(varargin)
% function mats = MATERIALS(varargin)

mats.MAT = cell(0,1);
for k=1:nargin
    if isa(varargin{k},'cell')
        mats_k = MATERIALS(varargin{k}{:}) ;
        mats.MAT = [mats.MAT , mats_k.MAT];
    elseif isa(varargin{k},'ELEMENT')
        mats.MAT = [mats.MAT {getmaterial(varargin{k})}];
    elseif isa(varargin{k},'MATERIAL')
        mats.MAT = [mats.MAT {varargin{k}}];
    elseif isa(varargin{k},'MATERIALS')
        mats.MAT = [mats.MAT varargin{k}.MAT];
    end
end
mats.n = length(mats.MAT);

mats = class(mats,'MATERIALS');

mats = setnumber(mats);
mats = unique(mats);
