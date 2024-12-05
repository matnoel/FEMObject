function u = solve_alphaupdate(A,b,W,Wtilde)

W = gathervectors(W);
W.alpha=1;
if nargin==3
Wtilde = W;
else
Wtilde = gathervectors(Wtilde);
Wtilde.alpha=1;    
end

WAW = timesblock(Wtilde'*A*W);
Wb =  timesblock(Wtilde'*b);
W.alpha = full((WAW\Wb)');
u=splitvectors(W);

