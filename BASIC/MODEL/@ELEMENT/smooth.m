function se = smooth(elem,node,se,gauss,varargin)
% function se = smooth(elem,node,se,gauss,varargin)

nbelem = getnbelem(elem);

if ischarin('ddl',varargin)
    ddl = getcharin('ddl',varargin,getddlgaussdual(elem));
    nbddlpergauss = length(ddl);
elseif ischarin('type',varargin)
    type = getcharin('type',varargin,'ddlgaussdual');
    nbddlpergauss = getnbddlpernode(elem,type);
    ddl = getddl(elem,type);
else
    nbddlpergauss = getcharin('nbddlpergauss',varargin,getnbddlpergauss(elem));
    ddl = getddlgaussdual(elem);
end
% nbddlpergauss = getnbddlpergauss(elem);
% nbddlpergauss = size(se,1);
n = size(se,2);
xnode = getcoord(node,getconnec(elem)');

% [Ns,detJ] = calc_Ngeom(elem,xnode,gauss.coord,'type','ddlgaussdual');
[Ns,detJ] = calc_Ngeom(elem,xnode,gauss.coord,'nbddlpernode',nbddlpergauss);
smd = sum(gauss.w*detJ*Ns'*se,4);

if getcharin('display',varargin,0)
    fprintf('\n smoothing\n')
end

M = MODEL(getmode(getsyscoord(node)));
M = addnode(M,node);
M = addelem(M,elem);
M = createddlnode(M,ddl);
masse = calc_massgeom(M);

elem2 = actualise_ddl(elem,ddl);

smd = assemble_vectorelem(M,{smd});
smd = masse\smd;
rep = getpos(node,unique(getconnec(elem)));
smd = reshape(full(smd),[nbddlpergauss,length(rep),n]);

smd = permute(smd,[1,3,2]);
se = zerosND(nbddlpergauss,n,getnbnode(node));
se(:,:,rep) = smd;
