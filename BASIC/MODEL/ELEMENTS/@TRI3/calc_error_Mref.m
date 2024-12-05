function [norms,fieldtype,fieldstorage,fieldddl] = calc_error_Mref(group,node,Mref,uref,u,varargin)
% function [norms,fieldtype,fieldstorage,fieldddl] = calc_error_Mref(group,node,Mref,uref,u,varargin)

fieldtype = 'ddlgaussdual';
fieldstorage = 'center';
fieldddl = getddlgaussdual(group);
PCref = getPC(uref);
lspcval = getvalue(Mref.ls{1});

norms{1} = zerosND(1,1,getnbelem(group));
norms{2} = zerosND(1,1,getnbelem(group));

for j=1:getnbelem(group)
    
    % fprintf('elem %d \n',j)
    elem = getelem(group,j);
    %ls=getlevelset(Mref.ls,getlsnumber(elem));
    n = getcharin('intorder',varargin,1);
    ng = getcharin('intordersto',varargin,getorder(PCref)+2);
    
    connec = getconnec(elem);
    xnode = node(connec);
    xnode = getcoord(xnode,connec');
    stomesh = getparam(elem,'stomesh');
    gausssto = calc_gausspoints(stomesh{1},ng);
    norms_elem{1} = zeros(gausssto.nbgauss,1);
    norms_elem{2} = zeros(gausssto.nbgauss,1);
    xg = transfer(RANDVARS(stomesh{1}),RANDVARS(PCref),gausssto.coord);
    u1 = randomlimiteval(u,xg,Mref); %problème mémoire
    uref1 = randomeval(uref,xg,Mref);
    lspce = lspcval(connec);
    lsk = randomeval(lspce,gausssto.coord,RANDVARS(stomesh));
    detJ = abs(calc_detJ(elem,xnode,gausssto.coord(1,:)));
    for e=1:gausssto.nbgauss
        if isa(u1,'MULTIMATRIX')
            ue = localize(elem,u1{e});
        elseif isa(u1,'double')
            ue = localize(elem,u1);
        end
        if isa(uref1,'MULTIMATRIX')
            ueref = localize(elem,uref1{e});
        elseif isa(u1,'double')
            ueref = localize(elem,uref1);
        end
        if isa(lsk,'double')
            [subgaussin,subgaussout] = lssubgauss_oneelem(lsk,n,2);
        else
            [subgaussin,subgaussout] = lssubgauss_oneelem(lsk{e},n,2);
        end
        if isempty(subgaussin.w)
        else
            [norms_elem{1}(e) norms_elem{2}(e)] = suberrorelem(elem,xnode,subgaussin,ue,ueref,detJ,varargin{:});
        end
    end
    norms{1}(1,1,j) = sum(norms_elem{1}.*gausssto.w);
    norms{2}(1,1,j) = sum(norms_elem{2}.*gausssto.w);
end
