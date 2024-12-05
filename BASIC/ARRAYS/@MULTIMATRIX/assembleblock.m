function w=assembleblock(u)

if isa(u.value,'cell')
    s=u.s;
    sm=u.sm;
    
    val = [];
    I = [];
    J = [];
    for k=1:sm(1)
        for l=1:sm(2)
            [i,j,v] = find(u.value{k,l});
            i = i+(k-1)*s(1);
            j = j+(l-1)*s(2);
            I = [I;i];
            J = [J;j];
            val = [val;v];
        end
    end
    w=sparse(I,J,val,s(1)*sm(1),s(2)*sm(2));
    
else
    s=u.s;
    sm=u.sm;
    
    ip=[];
    jp=[];
    vp=[];
    
    for i=1:sm(1)
        for j=1:sm(2)
            %    if isa(u.value,'double')
            temp = u.value(:,(j-1)*sm(1)+i);
            %    else
            %temp = u.value{i,j};
            %    end
            [ii,jj,v] = find(reshape(temp,s));
            ip=[ip;(i-1)*s(1)+ii];
            jp=[jp;(j-1)*s(2)+jj];
            vp=[vp;v];
        end
    end
    
    w=sparse(ip,jp,vp,s(1)*sm(1),s(2)*sm(2));
end