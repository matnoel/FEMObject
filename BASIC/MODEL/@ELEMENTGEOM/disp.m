function disp(u)
% function disp(u)

disp('Connectivity = [numelem , Inf , connec]')
disp(' ')
v = [u.numelem,Inf*ones(u.nbelem,1),u.connec] ;
disp(double(v))
disp([class(u) ' -> ELEMENTGEOM'])
v = struct('dim',getdim(u),'nbelem',getnbelem(u),'nbnode',getnbnode(u),'syscoordlocal',getsyscoordlocal(u),'numelem',getnumber(u),'connec',getconnec(u));
disp(v)
v = struct('option',getoption(u));
if ~isempty(v.option) && ~strcmp(v.option,' ')
    disp(v)
end
if getnbddl(u)>0
    clear v
    v.ddlnode = getddlname(getddlnode(u));
    v.ddlgauss = getddlname(getddlgauss(u));
    v.ddlnodedual = getddlname(getddlnodedual(u));
    v.ddlgaussdual = getddlname(getddlgaussdual(u));
    disp(v)
end

%if ~isempty(getlsnumber(u))
clear v
v.lsnumber = getlsnumber(u);
v.lsenrich = getlsenrich(u);
v.lstype = getlstype(u);
v.lsnature = getlsnature(u);
disp(v)
%end

v=getparam(u);
if length(fieldnames(v))
    disp('Additional parameters :')
    disp(v)
end


v=getmaterial(u);
if isa(v,'MATERIAL')
    disp(v);
else
    fprintf('  NO MATERIAL DEFINED\n');
end
