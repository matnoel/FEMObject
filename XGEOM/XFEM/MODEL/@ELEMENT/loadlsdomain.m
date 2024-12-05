function fe = loadlsdomain(elem,node,ls,compo,fun,varargin)
% function fe = loadlsdomain(elem,node,ls,compo,fun,varargin)
% compo : nom des composantes de la force
% fun : pointeur sur une fonction donnant la valeur des composantes en
% colonnes, peut N colonnes correspondant N cas de chargement
% appel de fun(x,varargin{:}) ou x sont le tableau de coordonnees de P points ou on veut evaluer la fonction
% size(x,1)=P size(x,2)=dimension de l'espace de definition des points
% la sortie est un 3D array de taille length(compo)*N*P  

n=getcharin('intorder',varargin,'lsload'); 

xnode = getcoord(node,getconnec(elem)');
ftest=fun(double(xnode(1,:,1)),varargin{:});

switch getlstype(elem)
    case 'out'
        fe = zerosND(elem.nbddl,size(ftest,2),getnbelem(elem));
    case {'in','indomain'}
        fe = load(elem,node,compo,fun,varargin{:});
    case 'cut'
        fe=zerosND(elem.nbddl,size(ftest,2),getnbelem(elem));
for e=1:getnbelem(elem)
    eleme = getelem(elem,e);
    connece = getconnec(eleme);
    xnodee = getcoord(node,getconnec(eleme)');
    lse = getvalue(ls);lse = lse(connece);
    if lsisin(eleme,lse)
        gauss = calc_gauss(eleme,n);
        fe(:,:,e) = fe(:,:,e) + integrate_with_gauss(eleme,xnodee,gauss,@eval_fe,compo,fun,varargin{:});
    elseif lsiscut(eleme,lse)
        [gaussin,gaussout] = calc_lssubgauss(eleme,lse,n);
        fe(:,:,e) = fe(:,:,e) + integrate_with_gauss(eleme,xnodee,gaussin,@eval_fe,compo,fun,varargin{:});
    end
end
end




function fe = eval_fe(xi,elem,xnode,compo,fun,varargin)
if getnbelem(elem)>0
[repddl,repddlcompo] = findddl(getddlnodedual(elem),compo);
N=calc_N(elem,xnode,xi);
x=calc_x(elem,xnode,xi);
f=fun(x,varargin{:});
f=MYDOUBLEND(f);
fe = N(repddl,:)'*f(repddlcompo,:);
else
ftest=fun(double(xnode(1,:)),varargin{:}) ;   
fe=zerosND(elem.nbddl,size(ftest,2),getnbelem(elem));
end

return

