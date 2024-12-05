function [Nls,detJ,x]=calc_Nls(elem,xnode,xgauss,ls,choix)


if nargout>1
detJ = calc_detJ(elem,xnode,xgauss);
end
if nargout>2
x=calc_x(elem,xnode,xgauss);
end


N = calc_N(elem,xnode,xgauss,'nbddlpernode',getindim(elem));

if getlsenrich(elem)==0
Nls=N;
else
    
lsxnode = getvalue(ls,getconnec(elem)');    
lsxnode = reshape(lsxnode,[size(lsxnode,1),1,size(lsxnode,2)]);
lsxnode = MYDOUBLEND(lsxnode);

Nuni = getN(elem,xgauss);

 signphix = sign(Nuni*lsxnode);
 signedlsxnode = signphix.*lsxnode;
 typeenrich=getenrichtype(ls);
 
switch typeenrich
    case 1
 lsx = Nuni*(abs(lsxnode)-signedlsxnode);
    case 2
 lsx = Nuni*signedlsxnode;   
    case 3
 global beta
            if isempty(beta); beta=0;  
            end
            connecfictitious = getparam(elem,'connecfictitious');
            if any(connecfictitious(:));
            factnode = ~connecfictitious;
            factnode = MYDOUBLEND(permute(factnode,[2,3,1]));  
            else
            factnode = 1;    
            end

lsx = (beta-signedlsxnode).*factnode;

    case 4
        global funenrich
 x = calc_x(elem,xnode,xgauss);
    lsx = funenrich(x);
end
        
%if strcmp(choix,'out')
%    lsx = Nuni*(abs(lsxnode)-lsxnode);
%elseif strcmp(choix,'in')
%    lsx = Nuni*(abs(lsxnode)+lsxnode);  
%else 
%    error('indiquer si c''est in ou out la levelset')
%end
    
Nls = [N , lsx*N];
rep = zeros(getindim(elem),getnbnode(elem),2);
rep(:)=1:numel(rep);
rep = permute(rep,[1,3,2]);
Nls(:,:) = Nls(:,rep(:)) ; 
end

