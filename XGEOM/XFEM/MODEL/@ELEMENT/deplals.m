function [se,fieldtype,fieldstorage,fieldddl] = deplals(elem,node,ls,q,varargin)
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

ls=getlevelset(ls,getlsnumber(elem));

if length(ls)>1
   error('pas programme pour une multilevelset') 
end

xnode = getcoord(node,getconnec(elem)');

qe = localize(elem,q);

switch getlstype(elem)
    case {'in','out'}
    if ~isempty(getmaterial(elem))
    N = calc_N(elem,xnode,gauss.coord);
    se = N*qe;
    else
    se=zerosND(elem.nbddlpergauss,size(qe,2),getnbelem(elem),numel(ls));  
    end
    case 'cut'
    N = calc_N(elem,xnode,gauss.coord);
    se = N*qe;
    
    if ~ischarin('calcout',varargin)
    N = calc_N(elem,xnode,gauss.coord,'nbddlpernode',1);
    
    lse = localize(elem,double(ls),'scalar');

    lse = N*lse;
    fact = sign(lse)<=0;
    se = se.*fact;
    end
        
end
if strcmp(getlstype(elem),'out')
se(:)=0;      
end
if ischarin('smooth',varargin)
se = smooth(elem,node,se,gauss);
end














