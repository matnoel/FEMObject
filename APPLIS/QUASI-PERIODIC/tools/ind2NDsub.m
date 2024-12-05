function sub = ind2NDsub(N,ind)

assert(max(ind)<=prod(N),'Size mismatch');

if issparse(ind)
    sparsity = true ;
    ind = full(ind) ;
else
    sparsity = false ; 
end

if isfloat(ind)
    ind = int64(ind) ;
end

d = numel(N) ;
sub = zeros(length(ind),d) ;

for i = 1:d
    c = kron( (1:N(i))' , ones(prod(N(1:i-1)),1) ) ;
    c = repmat(c,prod(N(i+1:end)),1) ;
    if sparsity
        sub(:,i) = sparse(c(ind)) ;
    else
        sub(:,i) = c(ind) ;
    end
end
%TODO : avoid creating full 'c'.

end