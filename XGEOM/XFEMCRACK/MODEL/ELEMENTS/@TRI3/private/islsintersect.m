function [rep,l2onl1] = islsintersect(l1,l2)

    segcon = [1,2,3;2,3,1]';
    seglsval = l1(segcon);   
    segcut = sign(prod(seglsval,2));
    repsegcut = find(segcut==-1);
    repsomcut = find(sign(l1)==0);
    if length(repsegcut)==1 && length(repsomcut)==1
    % levelset coupe 1 segment et passe par un sommet
    n = segcon(repsegcut,:);
    xi = l1(n(1))./(l1(n(1))-l1(n(2)));
    l21 = xi*l2(n(1)) + (1-xi)*l2(n(2));
    l22 = l2(repsomcut);
    rep = (sign(l21*l22)<=0);
    l2onl1 = [l21;l22];
    elseif length(repsegcut)==2 && isempty(repsomcut)
    % levelset coupe 2 segements
    n1 = segcon(repsegcut(1),:);
    n2 = segcon(repsegcut(2),:);    
    xi1 = l1(n1(1))./(l1(n1(1))-l1(n1(2)));
    xi2 = l1(n2(1))./(l1(n2(1))-l1(n2(2)));
    l21 = xi1.*l2(n1(1)) + (1-xi1).*l2(n1(2));
    l22 = xi2.*l2(n2(1)) + (1-xi2).*l2(n2(2));
    rep = (sign(l21*l22)<=0);
    l2onl1 = [l21;l22];
    
    else
       error('pas prevu'); 
    end