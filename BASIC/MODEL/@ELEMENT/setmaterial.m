function elem = setmaterial(elem,mat)
% function elem = setmaterial(elem,mat)

if isa(mat,'MATERIALS')
    error('on ne peut associer un MATERIALS a un ELEMENT -> rentrer un MATERIAL')
end

if ~isempty(mat)
    % if isa(mat,'MATERIAL') & isempty(getnumber(mat))
    %     error('on doit rentrer un material numerote')
    % end
end

elem.material = mat;