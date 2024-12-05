function tau = mylocelem(tau,elem)

if numel(tau)==1
    tau = MYDOUBLEND(tau);
else 
    tau = tau(getnumber(elem),:);
    tau = reshape(full(tau)',size(tau,2),1,getnbelem(elem));
    tau = MYDOUBLEND(tau);
end
