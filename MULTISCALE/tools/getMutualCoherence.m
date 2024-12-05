function mutual_coherence = getMutualCoherence(H,varargin)
% function mutual_coherence = getMutualCoherence(H,varargin)
% Computes the mutual coherence of matrix H
% mu(H) = max_{1<=i,j<=size(H,2),i~=j} |h_i'*h_j|/(||h_i||_2*||h_j||_2),
% where h_i is the ith column of matrix H

X = H/diag(sqrt(diag(H'*H)));
M = abs(X'*X);
M = M-diag(diag(M));
mutual_coherence = max(max(M));

end
