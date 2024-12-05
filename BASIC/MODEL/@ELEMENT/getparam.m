function param = getparam(elem,p)
% function param = getparam(elem,p)
% obtenir les parametres des elements p

param = elem.param ;
if nargin==2 && isa(p,'char')
    param = getfield(param,p);
elseif nargin==2
    
    f = fieldnames(param);
    param = struct2cell(param);
    
    param = [f(:),param(:)]';
    
    for i=1:length(f)
        if (isa(param{2,i},'MYDOUBLEND') || isa(param{2,i},'double')) &&  size(param{2,i},3)==getnbelem(elem)
            param{2,i} = getsubmatrix(param{2,i},p,3);
        elseif isa(param{2,i},'cell') && length(param{2,i}) == getnbelem(elem)
            param{2,i} = param{2,i}(p);
        elseif (isa(param{2,i},'double') || isa(param{2,i},'logical')) && (size(param{2,i},1) == getnbelem(elem))
            param{2,i} = param{2,i}(p,:);
        end
    end
    for i=1:length(f)
        if isa(param{2,i},'cell')
            param{2,i} = {param{2,i}};
        end
    end
    
    param = struct(param{:});
end