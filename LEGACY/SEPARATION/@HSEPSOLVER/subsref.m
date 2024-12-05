function s = subsref(A,s)
% function s = subsref(A,s)

if strcmp(s(1).type,'()')
    param = getparam(A);
    node  = s(1).subs{1};
    param.node           = node;
    param.maxorder       = param.maxorder       (node);
    param.maxiter        = param.maxiter        (node);
    param.itercrit       = param.itercrit       (node);
    param.tol            = param.tol            (node);
    param.display        = param.display        (node);
    param.dyntol         = param.dyntol         (node);
    param.dynitercrit    = param.dynitercrit    (node);
    param.residual       = param.residual       (node);
    param.errorindicator = param.errorindicator {node};
    param.alphaupdate    = param.alphaupdate    (node);
    param.update         = param.update         (node);
    param.updatetucker   = param.updatetucker   (node);
    
    c=getconnect(param.tree);
    N=node;
    while c(N)~=0
        param.depth      = param.depth +1 ;
        N=c(N);
    end
    s = cell2struct(struct2cell(param),fieldnames(param),1);
end
