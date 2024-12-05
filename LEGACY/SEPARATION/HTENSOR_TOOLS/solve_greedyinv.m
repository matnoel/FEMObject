function P = solve_greedyinv(A,rank1,errcrit,opts)
% function P = solve_greedyinv(A,rank1,errcrit,opts)

P = rank1(A,[]);

err = zeros(opts.maxrank,1);

for i = 1:opts.maxrank
    P0 = P;
    W = rank1(P);
    P = P+W;
    err(i) = errcrit(P,P0,W);
end

end

