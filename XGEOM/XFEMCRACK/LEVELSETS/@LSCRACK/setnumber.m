function c = setnumber(c,number)
% function c = setnumber(c,number)

if ~isa(number,'double')
    error('le numero doit etre un double')
end
c.number = number;

