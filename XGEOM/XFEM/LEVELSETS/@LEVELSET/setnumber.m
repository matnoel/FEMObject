function ls = setnumber(ls,number)
% function ls = setnumber(ls,number)

if ~isa(number,'double')
    error('le numero doit etre un double')
end
ls.number = number;