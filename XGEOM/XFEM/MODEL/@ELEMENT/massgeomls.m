function ke = massgeomls(elem,node,ls,varargin)

ls=getlevelset(ls,getlsnumber(elem));


n=getcharin('intorder',varargin,'mass');
typeenrich=getlsenrichtype(elem);
tic
switch getlstype(elem)
    case 'out'
   me = zerosND(getnbnode(elem),getnbnode(elem),getnbelem(elem),numel(ls));
    case {'in','indomain'}
   me = massgeom(elem,node);     
   me = repmat(me,[1,1,1,numel(ls)]);

    case 'cut'
   me = zerosND(getnbnode(elem),getnbnode(elem),getnbelem(elem),numel(ls));
xnode = getcoord(node,getconnec(elem)');

for i=1:numel(ls)
[elemin,elemcut,elemout,repin,repcut,repout,xnodein,xnodecut,xnodeout] = lssplitelem(elem,ls{i},node);

ke(:,:,repin,i) =ke(:,:,repin,i)+integrate(elemin,xnodein,n,@eval_me);

[elemcutin,elemcutout,nodeplus,xnodecutin,xnodecutout]=lsdivideelem(elemcut,ls{i},node);


ke(:,:,repcut,i)=ke(:,:,repcut,i)+...
    subintegrate(elemcut,xnodecut,elemcutin,xnodecutin,n,[getnbnode(elem),getnbnode(elem)],@eval_me);

end

fprintf('Elapsed time is %.3f seconds.',toc)


end


function me = eval_me(xi,elem,xnode)

if nargin>3 & getnbelem(elem)>0 & ~strcmpi(getlstype(elem),'out')
N=calc_Ngeom(elem,xnode,xi,varargin{:});
me=N'*N;
else
me=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));   
end
return

