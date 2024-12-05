function mat = setnumber(mat,number)
% function mat = setnumber(mat,number)

if ~isa(number,'double')
    error('The number must be a double.')
end
mat.number = number;
