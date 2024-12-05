function cellConn = calcMesoConnectivity(Assembler,direction)
% cellConn = calcMesoConnectivity(Assembler,direction)

cellNum = getCellNum(Assembler) ;
N = cellNum(direction) ;
orthDir = 1+mod(direction+2,2) ; % 2 if direction=1, 1 if direction=2

cellConn = sparse(1:N-1,2:N,ones(N-1,1),N,N) ;

% If periodic boundary condition, connect final cell to first cell
if isPeriodic(Assembler,direction)
    cellConn(N,1) = 1 ;
end

if getOrder(Assembler)==2 % then assemble global connectivity matrix
    factors = {eye(cellNum(orthDir)),cellConn} ;
    cellConn = kron(factors{direction},factors{orthDir}) ;
end
end