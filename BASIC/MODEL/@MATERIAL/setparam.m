function mat = setparam(mat,i,value)
% function mat = setparam(mat,i,value)

value = full(value);
if isa(i,'double')
    mat.param{i,2} = value;
elseif isa(i,'char') || isa(i,'cell')
    [rep,pos] = ischarin(i,mat.param(:,1));
    mat.param{pos,2} = value;
end

