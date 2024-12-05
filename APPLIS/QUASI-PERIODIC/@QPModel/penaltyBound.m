function penalty = penaltyBound(model,bounds,useSW,useAW,isPeriodic)
% penalty = penaltyBound(model,bounds,useSW,useAW,isPeriodic)

if nargin < 5
    if nargin < 4
        useAW = false ;
        if nargin < 3
            useSW = false ;
        end
    end 
    assert(~(useAW||useSW),['Cannot apply required weights: ',...
        'Periodicity missing.']) ;
    isPeriodic = [false false] ;
elseif numel(isPeriodic)==1
    isPeriodic = [isPeriodic isPeriodic] ;
end

traceCst = getTraceConstantL2(model) ;
if useSW || useAW
    cellConn = mesoConnectivity(2,getCellNum(model),1,isPeriodic(1));
    cellConn = cellConn | mesoConnectivity(2,getCellNum(model),2, ...
                                           isPeriodic(2)) ;
else
    cellConn = [] ;
end

penalty = penaltyBound(bounds,traceCst,getCellSize(model),...
    useSW,useAW,cellConn) ;
ifprint(model.verbose,'Infimum penalty %g\n',penalty)
end