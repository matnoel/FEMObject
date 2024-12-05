function R = project(v,f)
% function R = project(v,f)

R = v.MYDOUBLEND';

if nargin==2
    R = R*f;
end
