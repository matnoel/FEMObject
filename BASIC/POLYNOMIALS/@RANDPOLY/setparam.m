function h = setparam(h,field,val)
% function h = setparam(h,field,val)

if nargin==2
    h.param = field;
else
    h.param= setfield(h.param,field,val);  
end


