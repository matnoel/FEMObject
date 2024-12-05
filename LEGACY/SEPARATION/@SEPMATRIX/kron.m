function A = kron(Asep)

if Asep.dim==2
    
    A = Asep.alpha(1)*kron(Asep.F{1,1},Asep.F{1,2});
    for k=2:Asep.m
    A = A + Asep.alpha(k)*kron(Asep.F{k,1},Asep.F{k,2});    
    end
    
else
    
    error('not implemented')    
end