function u=ctranspose(u)

u.MYDOUBLE = u.MYDOUBLE';

if ndims(u)>2
    error('ctranspose non defini')
end

if isa(u,'FEVECTOR')
fe0 = u.fedims ;
fe1 = setdiff(1:2,fe0);
u.fedims = fe1 ; 
u.repbloque{fe0}=':';
u.repbloque{fe1}=u.ddlbloque;
u.repfree{fe0}=':';
u.repfree{fe1}=u.ddlfree;
end