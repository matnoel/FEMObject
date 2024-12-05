function [ T ] = orth( w,varargin )
% function [ T ] = orth( w )
% function [ T ] = orth( w,dim )
% Orthogonalize the TSEPMATRIX w

if nargin==1
    dim=1:w.dim;
else
    dim=varargin{1};
end

if size(w.F{1}{1},2)==1
    gw=gathervectors(w);
else
    gw=gathervectors(w');
end
for d=dim
    if(size(gw.F{d},2)~=1)
        qd=orth(full(gw.F{d}));
        rd=qd'*gw.F{d};
        gw.F{d}=qd;
        gw.alpha=ttm(gw.alpha,rd,d);
    end
end
gw.m=size(gw.alpha);
if size(w.F{1}{1},2)==1
    T=splitvectors(gw);
else
    T=splitvectors(gw)';
end

if any(T.m ~= size(T.alpha))
    sz=size(T.alpha);
    dim=find(T.m~=size(T.alpha));
    for d=dim
        n=T.m(d)/sz(d);
        f=cell(sz(d),1);
        for m=1:sz(d)
            f{m}=[T.F{d}{n*(m-1)+1:n*m}];
        end
        T.F{d}=f;
        T.m(d)=sz(d);
    end
end



end

