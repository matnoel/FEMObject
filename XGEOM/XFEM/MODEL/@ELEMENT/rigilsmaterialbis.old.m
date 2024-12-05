function ke = rigilsmaterial(elem,node,ls,varargin)

if isempty(getlsnumber(elem))
    ke = rigi(elem,node,varargin{:});
    return
end

ls=getlevelset(ls,getlsnumber(elem));
if ~isempty(getmaterial(ls))
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end

n=getcharin('intorder',varargin,'rigils');
typeenrich=getenrichtype(ls);

switch getlstype(elem)
case {'in','indomain','out'}
    switch getlstype(elem)
    case {'in'}
        elem = setmaterial(elem,matin);
    case {'out'}
        elem = setmaterial(elem,matout);
    end

    if ~isempty(getmaterial(elem))
        if typeenrich<=1
            ke = rigi(elem,node);
            ke=repmat(ke,[1,1,1,numel(ls)]);
        elseif typeenrich>1
            ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem),numel(ls));    
            xnode = getcoord(node,getconnec(elem)');
            for i=1:numel(ls)
                ke(:,:,:,i) = integrate(elem,xnode,n,@eval_ke,getmaterial(elem),ls{i},getlstype(elem));
            end
        end

    else
        ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem),numel(ls));    
    end

case 'cut'
    ke=zerosND(elem.nbddl,elem.nbddl,getnbelem(elem),numel(ls));
    xnode = getcoord(node,getconnec(elem)');

    for i=1:numel(ls)
        [elemin,elemcut,elemout,repin,repcut,repout,xnodein,xnodecut,xnodeout] = lssplitelem(elem,ls{i},node);

        ke(:,:,repin,i) =ke(:,:,repin,i)+integrate(elemin,xnodein,n,@eval_ke,matin,ls{i},'in');
        ke(:,:,repout,i)=ke(:,:,repout,i)+integrate(elemout,xnodeout,n,@eval_ke,matout,ls{i},'out');


        [elemcutin,elemcutout,nodeplus,xnodecutin,xnodecutout]=lsdivideelem(elemcut,ls{i},node);


        ke(:,:,repcut,i)=ke(:,:,repcut,i)+...
            subintegrate(elemcut,xnodecut,elemcutin,xnodecutin,n,[elem.nbddl,elem.nbddl],@eval_ke,matin,ls{i},'in');


        ke(:,:,repcut,i)=ke(:,:,repcut,i)+...
            subintegrate(elemcut,xnodecut,elemcutout,xnodecutout,n,[elem.nbddl,elem.nbddl],@eval_ke,matout,ls{i},'out');

    end

    fprintf('Elapsed time is %.3f seconds.',toc)


end


function ke = eval_ke(xi,elem,xnode,mat,ls,choix)

if nargin>3 & ~isempty(mat) & getnbelem(elem)>0 & ~strcmpi(getlstype(elem),'out')
    D=calc_opmat(mat,elem,xnode,xi);
    Bls=calc_Bls(elem,xnode,xi,ls,choix);

    ke = Bls'*D*Bls ;
else
    ke=zerosND(getnbddl(elem),getnbddl(elem),getnbelem(elem));
end


return

