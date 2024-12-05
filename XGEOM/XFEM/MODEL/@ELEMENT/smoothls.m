function se = smoothls(elem,node,se,gauss,ls,varargin)

nbddlpergauss = getnbddlpergauss(elem);
nbelem = getnbelem(elem);
n=size(se,2);
xnode = getcoord(node,getconnec(elem)');

[Ns,detJ] = calc_Ngeom(elem,xnode,gauss.coord,'type','ddlgaussdual');
smd = sum(gauss.w*detJ*Ns'*se,4);

if getcharin('display',varargin,0)
fprintf('\n smoothing\n')
end

M = MODEL(getmode(getsyscoord(node)));
M = addnode(M,node);
M = addelem(M,elem);
M = createddlnode(M,getddlgaussdual(elem));
masse=calc_massgeomls(M);


elem2=actualise_ddl(elem,getddlgaussdual(elem));  


smd = assemble_vectorelem(M,{smd});
smd = masse\smd;
rep = getpos(node,unique(getconnec(elem)));
smd = reshape(full(smd),[nbddlpergauss,length(rep),n]);

smd = permute(smd,[1,3,2]);
se=zerosND(nbddlpergauss,n,getnbnode(node));
se(:,:,rep)=smd;
