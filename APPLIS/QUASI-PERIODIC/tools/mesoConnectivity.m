function conn = mesoConnectivity(order,cellNum,direction,isPeriodic)
%  conn = mesoConnectivity(order,cellNum,direction,periodicity)

    N = cellNum(direction) ;

    % If periodic boundary condition, connect final cell to first cell
    if isPeriodic
        sub1 = (1:N)' ;
        sub2 = [2:N 1]' ;
    else
        sub1 = (1:N-1)' ;
        sub2 = (2:N)' ;
    end
    conn = sparse(sub1,sub2,ones(numel(sub1),1),N,N) ;

    if order==2 % then assemble global connectivity matrix
        orthDir = 1+mod(direction+2,2) ; % 2 if direction=1, 1 if direction=2
        factors = {eye(cellNum(orthDir)),conn} ;
        conn = kron(factors{direction},factors{orthDir}) ;
    end

end