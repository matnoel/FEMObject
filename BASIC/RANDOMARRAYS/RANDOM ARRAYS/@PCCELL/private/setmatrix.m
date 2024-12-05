function B = setmatrix(A,n1,n2)

for k=1:size(A,2)
B{k} = reshape(A(:,k),n1,n2);
end

return
