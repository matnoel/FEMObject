function var=addcharin(c,var)
% function var=addcharin(c,var)

if ~isa(c,'char')
    error('argument doit etre char')
end
if ~ischarin(c,var)
    var=[var , {c}]; 
end
