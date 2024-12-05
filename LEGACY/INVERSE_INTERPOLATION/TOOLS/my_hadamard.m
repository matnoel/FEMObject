function H = my_hadamard(n,m,H)

H0=[1 1;1 -1];
s0=[2 2];

if nargin==2
    H=H0;
    
    
end

s=size(H);

if (s(1)==n) && (s(2)==m)
    return
end

if (s0(1)*s(1)<=n) && (s0(2)*s(2)<=m)
    H=my_hadamard(n,m,kron(H0,H));
    return
end


if s(1)<n
    ind1=min((n-s(1)),s(1));
    H21=H0(2,1) * H(1:ind1,:);
else H21=[];
end

if s(2)<m
    ind2=min((m-s(2)),s(2));
    H12=H0(1,2) * H(:,1:ind2);
else H12=[];
end

if (s(1)<n) && (s(2)<m)
    H22=H0(2,2) * H(1:ind1,1:ind2);
else H22=[];
end

H=my_hadamard(n,m, [H H12 ;H21 H22] );


