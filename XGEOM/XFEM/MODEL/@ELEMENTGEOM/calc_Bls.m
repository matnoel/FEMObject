function [Bls,detJ]=calc_Bls(elem,xnode,xgauss,ls,choix)


Nuni = getN(elem,xgauss);
[DN,detJ] = calc_DN(elem,xnode,xgauss);

B = calc_B(elem,xnode,xgauss);

if getlsenrich(elem)==0
Bls= B ;

else

     
lsxnode = getvalue(ls,getconnec(elem)');    


typeenrich=getenrichtype(ls);
switch typeenrich
    case 1
if strcmp(choix,'out') || strcmp(choix,'indomain')
    lsx = Nuni*(abs(lsxnode)-lsxnode);
    dlsx = DN*(abs(lsxnode)-lsxnode);
elseif strcmp(choix,'in') 
    lsx = Nuni*(abs(lsxnode)+lsxnode);  
    dlsx = DN*(abs(lsxnode)+lsxnode);
end

    case 2
if strcmp(choix,'out') | strcmp(choix,'indomain')
    lsx = Nuni*lsxnode;
    dlsx = DN*lsxnode;
elseif strcmp(choix,'in')
    lsx = -Nuni*lsxnode;  
    dlsx = -DN*lsxnode;
end        
    
    case 3
            global beta
            if isempty(beta)
                beta=0;
            end
            connecfictitious = getparam(elem,'connecfictitious');
            if any(connecfictitious(:));
            factnode = ~connecfictitious;
            factnode = MYDOUBLEND(permute(factnode,[2,3,1]));  
            else
            factnode = 1;    
            end
if strcmp(choix,'out') | strcmp(choix,'indomain')
  tildenode = (beta-lsxnode).*factnode; %  valeurs nodales de beta-|phi| 
  lsx = Nuni*tildenode;
  dlsx = DN*tildenode;     
 elseif strcmp(choix,'in')
  tildenode = (beta+lsxnode).*factnode; %  valeurs nodales de beta-|phi| 
  lsx = Nuni*tildenode;
  dlsx = DN*tildenode;  
end        

    case 4
        global funenrich
        global dfunenrich
    x = calc_x(elem,xnode,xgauss);
    lsx = funenrich(x);
    
    dlsx = dfunenrich(x);
    
end



switch getindim(elem)
    case 1
Dlsx=dlsx;
    case 2
Dlsx = zerosND([3,2,sizeND(dlsx)]);
Dlsx(1,1) = dlsx(1);
Dlsx(2,2) = dlsx(2);
Dlsx(3,1) = dlsx(2);
Dlsx(3,2) = dlsx(1);
    case 3
Dlsx = [dlsx(1),0,0;0,dlsx(2),0;0,0,dlsx(3);dlsx(2),dlsx(1),0;dlsx(3),0,dlsx(1);0,dlsx(3),dlsx(2)];            
end

N = calc_N(elem,xnode,xgauss,'nbddlpernode',getindim(elem));

Bls = [B ,lsx*B + Dlsx*N];
rep = zeros(getindim(elem),getnbnode(elem),2);
rep(:)=1:numel(rep);
rep = permute(rep,[1,3,2]);
Bls(:,:) = Bls(:,rep(:)) ; 


end

