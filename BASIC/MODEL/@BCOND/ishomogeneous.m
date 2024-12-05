function r=ishomogeneous(B,i)

if nargin==1
    r=ones(1,length(B));
    for i=1:length(B.BC)
        r(i)=B.BC{i}.ishomogeneous;
    end
else
    r=ishomogeneous(B);
    r=r(i);
end
