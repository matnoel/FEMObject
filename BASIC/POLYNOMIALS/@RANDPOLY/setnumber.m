function rv = setnumber(rv,number)
% function rv = setnumber(rv,number)

if ~isa(number,'double')
    error('le numero doit etre un double')
end
rv.number = number;
