function A = TIMEMATRIX(R)

T = R.TIMEMODEL ;
s  = size(R.V);

A = double(R.V) * R.D * double(R.L);
A = TIMEMATRIX(A,s,T);