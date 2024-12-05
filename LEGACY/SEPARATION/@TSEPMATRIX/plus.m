function w = plus(u,v,varargin)

w = u;
w.m = u.m+v.m;
for i=1:u.dim
    w.F{i}=[u.F{i};v.F{i}];
end

if all(u.m==0)
    w=v;
elseif all(v.m==0)
    w=u;
else
    if w.dim==1
        w.alpha=tensor(zeros(w.m,1),w.m);
        w.alpha(1:u.m(1))=u.alpha;
        w.alpha((u.m(1)+1):end)=v.alpha;
    else
        w.alpha=tensor(zeros(w.m),w.m);
        switch w.dim
            case 2
                w.alpha(1:u.m(1),1:u.m(2))=u.alpha;
                w.alpha((u.m(1)+1):end,(u.m(2)+1):end)=v.alpha;
            case 3
                w.alpha(1:u.m(1),1:u.m(2),...
                    1:u.m(3))=u.alpha;
                w.alpha((u.m(1)+1):end,(u.m(2)+1):end,...
                    (u.m(3)+1):end)=v.alpha;
            case 4
                w.alpha(1:u.m(1),1:u.m(2),...
                    1:u.m(3),1:u.m(4))=u.alpha;
                w.alpha((u.m(1)+1):end,(u.m(2)+1):end,...
                    (u.m(3)+1):end,(u.m(4)+1):end)=v.alpha;
            case 5
                w.alpha(1:u.m(1),1:u.m(2),...
                    1:u.m(3),1:u.m(4),...
                    1:u.m(5))=u.alpha;
                w.alpha((u.m(1)+1):end,(u.m(2)+1):end,...
                    (u.m(3)+1):end,(u.m(4)+1):end,...
                    (u.m(5)+1):end)=v.alpha;
        end
    end
end
