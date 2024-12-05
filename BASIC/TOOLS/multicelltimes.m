function c = multicellmtimes(a,b)
%function c = multicellmtimes(a,b)

if ~isa(a,'cell') & ~isa(b,'cell')
    error('un argument doit etre une cellule')
elseif isa(a,'cell') & isa(b,'double')
    sa = size(a);
    sb = size(b);
    sda = size(a{1});
    if ~all(sa==sb) & ~all(sb==1)
        error('les dimensions ne correspondent pas')
    end

    if all(sb==1)
        c = a; 
        for k=1:numel(a)
            c{k}=sparsemtimes(a{k},b)  ;
        end
    else
        c=cell(sa(1),sa(2));
        for i=1:sa(1)
            for j=1:sa(2)
                c{i,j} = sparse(sda(1),sda(2)); 
            end
        end

        [J,K] = find(b);
        for n=1:length(J)   
            j = J(n);
            k = K(n);
            c{j,k} = a{j,k}*b(j,k);
        end

    end

elseif isa(b,'cell') & isa(a,'double')
    c = multicelltimes(b,a);
else
    error('pas prevu')
end




