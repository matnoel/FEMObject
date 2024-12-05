function [se,fieldtype,fieldstorage,fieldddl] = sigmals(elem,node,ls,q,varargin)


fieldtype = 'ddlgaussdual';
fieldddl = getddlgaussdual(elem);

if ischarin('node',varargin)
gauss.coord=permute(nodelocalcoord(elem),[4,2,3,1]);
gauss.nbgauss = getnbnode(elem);
fieldstorage = 'node';
elseif ischarin('smooth',varargin)
n=getcharin('intorder',varargin,'mass');
gauss=calc_gauss(elem,n);        
fieldstorage = 'node';
else
n=getcharin('intorder',varargin,'rigi');
gauss=calc_gauss(elem,n);        
fieldstorage = 'gauss';
end

if isempty(getlsnumber(elem))
    matin = getmaterial(elem);
    matout = [];
else
 ls=getlevelset(ls,getlsnumber(elem));

if length(ls)>1
   error('pas programme pour une multilevelset') 
end

if ~isempty(getmaterial(ls))
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end
end


xnode = getcoord(node,getconnec(elem)');

qe=localize(elem,q);

switch getlstype(elem)
    case {'in','indomain','out'}

switch getlstype(elem)
    case {'in'}
        elem = setmaterial(elem,matin);
    case {'out'}
        elem = setmaterial(elem,matout);
end

if ~isempty(getmaterial(elem))
se = sigma(getmaterial(elem),elem,xnode,gauss.coord,qe);
else
se=zerosND(elem.nbddlpergauss,size(qe,2),getnbelem(elem),numel(ls));    
end

    case {'cut','touchcut'}
switch getnature(ls)
    case 'domain'
elem = setmaterial(elem,matin);
se = sigma(getmaterial(elem),elem,xnode,gauss.coord,qe);
if ischarin('nocalcout',varargin)
N = calc_N(elem,xnode,gauss.coord,'nbddlpernode',1);  
lse = localize(elem,double(ls),'scalar');
lse = N*lse;
fact = sign(lse)<=0;
se = se.*fact;
end
    case 'material'
        % pas bien fait  
seout = eval_se(gauss.coord,elem,xnode,qe,getmaterial(elem),ls,'out');
sein = eval_se(gauss.coord,elem,xnode,qe,getmaterial(ls),ls,'in');
N = calc_N(elem,xnode,gauss.coord,'nbddlpernode',1);        
lse = localize(elem,double(ls),'scalar');
lse = N*lse;
factin = sign(lse)<=0;
factout = sign(lse)>0;
se = sein.*factin + seout.*factout;
end
end


if ischarin('smooth',varargin)
se = smoothls(elem,node,se,gauss,ls);
end




function se = eval_se(xi,elem,xnode,qe,mat,ls,choix)
if getnbelem(elem)>0
D=calc_opmat(mat,elem,xnode,xi);
Bls=calc_Bls(elem,xnode,xi,ls,choix);
se =D*Bls*qe ;
else
se=zerosND(elem.nbddlpergauss,size(qe,2),getnbelem(elem));
end

return


