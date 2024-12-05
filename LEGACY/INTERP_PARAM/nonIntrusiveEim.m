function [approx,R,lambda,Xim,indx] = nonIntrusiveEim(N,fun,snapshotsParamFun,Xi,optsEim)
% function [approx,R,lambda,Xim,indx] = nonIntrusiveEim(N,fun,snapshotsParamFun,Xi,optsEim)
    [~,~,indxi,indx] = eim(snapshotsParamFun,optsEim);
    lambda = snapshotsParamFun(:,indxi);
    Xim = Xi(:,indxi);
    nsnap = numel(indxi);
    R = zeros(N,nsnap);
    for k = 1:nsnap
        R(:,k) = fun(Xim(:,k));
    end
    lambda = lambda(indx,:) \ snapshotsParamFun(indx,:);
    approx = R*lambda;
end
