function a = getvalue(a)
% function a = getvalue(a)

if getnbdim(a)==0
    a = a.factor;
elseif getnbdim(a)==1
    a = a.phi{1}*a.factor;
else
    error('getvalue n''a pas de sens ')
end