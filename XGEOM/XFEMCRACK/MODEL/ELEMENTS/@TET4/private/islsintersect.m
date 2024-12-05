function [rep,l2onl1] = islsintersect(l1,l2)

   
    segcon = [1,2;2,3;1,3;1,4;2,4;3,4];
    seglsval = l1(segcon);   
    segcut = sign(prod(seglsval,2));
    repsegcut = find(segcut==-1);
    repsomcut = find(sign(l1)==0);
    
    if length(repsegcut)==3 && isempty(repsomcut)
        l2onl1 = zeros(1,3);
     for k=1:3
         n=segcon(repsegcut(k),:);
         xi = l1(n(1))./(l1(n(1))-l1(n(2)));
         l2onl1(k) = (1-xi)*l2(n(1)) + (xi)*l2(n(2));
     end
%     rep = any(sign(prod(l2onl1([1,2;1,3;2,3]),2))<=0);
     rep = any(l2onl1>0) && any(l2onl1<=0);
    elseif length(repsegcut)==4 && isempty(repsomcut)
         l2onl1 = zeros(1,4);
     for k=1:4
         n=segcon(repsegcut(k),:);
         xi = l1(n(1))./(l1(n(1))-l1(n(2)));
         l2onl1(k) = (1-xi)*l2(n(1)) + (xi)*l2(n(2));
     end
     
  %   rep = any(sign(prod(l2onl1([1,2;1,3;1,4;2,3;2,4;3,4]),2))<=0);
    rep = any(l2onl1>0) && any(l2onl1<=0);
    elseif length(repsegcut)==2 && length(repsomcut)==1
        l2onl1 = zeros(1,3);
        for k=1:2
         n=segcon(repsegcut(k),:);
         xi = l1(n(1))./(l1(n(1))-l1(n(2)));
         l2onl1(k) = (1-xi)*l2(n(1)) + (xi)*l2(n(2));
        end
        l2onl1(3) = l2(repsomcut);
        %rep = any(sign(prod(l2onl1([1,2;1,3;2,3]),2))<=0);
        rep = any(l2onl1>0) && any(l2onl1<=0);
    elseif length(repsegcut)==1 && length(repsomcut)==2
        l2onl1 = zeros(1,3);
        for k=1:1
         n=segcon(repsegcut(k),:);
         xi = l1(n(1))./(l1(n(1))-l1(n(2)));
         l2onl1(k) = (1-xi)*l2(n(1)) + (xi)*l2(n(2));
        end
        l2onl1(2:3) = l2(repsomcut);
        %rep = any(sign(prod(l2onl1([1,2;1,3;2,3]),2))<=0);
        rep = any(l2onl1>0) && any(l2onl1<=0);
    else
        keyboard
       error('pas prevu'); 
    end