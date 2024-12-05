function mat = MATERIAL(varargin)

if nargin==0
    mat.type = '';
    mat.param = cell(0,2);
    mat.number = [];
    mat = class(mat,'MATERIAL');
elseif isa(varargin{1},'MATERIAL')
    mat = varargin{1};
elseif isa(varargin{1},'ELEMENT')
    mat = getmaterial(varargin{1});
else
    mat.type = varargin{1};
    
    names = fieldnames(varargin{2});
    values = struct2cell(varargin{2});
    
    number = getcharin('number',varargin);
    
    mat.param = cell(length(names),2);
    mat.param(:,1) = names;
    mat.param(:,2) = values;
    mat.number = number;
    
    mat = class(mat,'MATERIAL');
end
