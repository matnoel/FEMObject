function w=mtimes(u,v)

if isa(u,'FEVECTOR') & isa(v,'FEMATRIX')
    w=FEVECTOR(MYDOUBLE(u)*v.MYDOUBLE,u.ddlbloque,2);
elseif isa(u,'FEMATRIX') & isa(v,'FEVECTOR')   
    w=FEVECTOR(u.MYDOUBLE*MYDOUBLE(v),u.ddlbloque,1);
elseif isa(u,'FEMATRIX') & isa(v,'FEMATRIX')
    w=FEMATRIX(u.MYDOUBLE*v.MYDOUBLE,u.ddlbloque);
elseif (isa(u,'FEMATRIX')|isa(u,'FEVECTOR')) & (isa(v,'MYDOUBLE')|isa(v,'double'))   
    w=u;
    w.MYDOUBLE = u.MYDOUBLE * MYDOUBLE(v) ;
elseif (isa(v,'FEMATRIX')|isa(v,'FEVECTOR')) & (isa(u,'MYDOUBLE')|isa(u,'double'))   
    w=v;
    w.MYDOUBLE = MYDOUBLE(u) * v.MYDOUBLE ;
end

