function T = progressive_tucker_als(X,allR,tol)
% function T = progressive_tucker_als(X,allR,tol)

if isa(X,'double')
    X=tensor(X);
end

for R=allR
    [T,fit]=tucker_als(X,R);
    ctol=1-fit; %erreur norm L2
    if ctol<tol
        break;
    end
end

if isa(T,'ttensor')
    T=TSEPMATRIX(T);
end

end
