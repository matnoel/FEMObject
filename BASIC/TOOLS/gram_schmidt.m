function W = gram_schmidt(W,M)
% function W = gram_schmidt(W,M)
% Gram schmidt orthogonalization of matrix W
% M is the matrix defining the inner product

for i=1:size(W,2)
    for j=1:i-1  
        W(:,i) = W(:,i)-W(:,j)*(W(:,j)'*M*W(:,i)); 
    end
    n = sqrt(W(:,i)'*M*W(:,i));
    W(:,i)=W(:,i)/n;
end


