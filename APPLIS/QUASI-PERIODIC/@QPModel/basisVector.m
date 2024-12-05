function v = basisVector(model,dir)
% v = basisVector(model,dir)
% Returns sign(dir)*e_abs(dir)

N = 2*getNbCellDoF(model) ;
v = distributeMicro(model,{(1:getCellNb(model))'},...
    sparse(abs(dir):2:N,1,sign(dir),N,1)) ;
end