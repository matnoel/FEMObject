function c = multicellmtimes(a,b)
%function c = multicellmtimes(a,b)

if ~isa(a,'cell') & ~isa(b,'cell')
    error('un argument doit etre une cellule')
elseif isa(a,'cell') & isa(b,'double')
    sa = size(a);
    sb = size(b);
    sda = size(a{1});

    if sa(2)~=sb(1) & ~all(sb==1)
        error('les dimensions ne correspondent pas')
    end

    if all(sb==1)
        c=cell(sa(1),sa(2));    
        for i=1:sa(1)
            for j=1:sa(2)
                c{i,j} = sparsemtimes(a{i,j},b); 
            end
        end 

    else


        c=cell(sa(1),sb(2));
        for i=1:sa(1)
            for j=1:sb(2)
                c{i,j} = sparse(sda(1),sda(2)); 
            end
        end

        [J,K] = find(b);
        for n=1:length(J)
            for i=1:sa(1);
                j = J(n);
                k = K(n);  
                c{i,k} = c{i,k} +  sparsemtimes(a{i,j},b(j,k));
            end
        end

    end

elseif isa(b,'cell') & isa(a,'double')
    c = multicellmtimes(b',a')';

else
    error('pas prevu')
end




