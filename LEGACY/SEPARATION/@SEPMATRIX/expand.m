function A = expand(u)


if u.dim==1 && u.m>0 && size(u.F{1,1},2)==1
    u = gathervectors(u);
    A = u.F{:,1}*diag(u.alpha);
elseif all(size(u)==1)
    A = reshape([u.F{:}],size(u.F,1),size(u.F,2));
    A = full(sum(u.alpha(:).*prod(A,2)));

elseif u.dim==1
    A = u.F{1,1}*u.alpha(1);
    for j=2:u.m
        A = A + u.F{j,1}*u.alpha(j);
    end
elseif u.dim==2
    mm=length(u.alpha);
    u = gathervectors(u);
    A = u.F{1,1}*spdiags(u.alpha(:),0,mm,mm)*u.F{1,2}';
elseif u.dim==3
    mm=length(u.alpha);
    u = gathervectors(u);
    s = size(u);
    A = zeros(s);
    for l=1:s(3)
        d=u.alpha.*u.F{1,3}(l,:);
        A(:,:,l) = u.F{1,1}*(spdiags(d(:),0,mm,mm))*u.F{1,2}';
    end
elseif u.dim==4
    mm=length(u.alpha);
    u = gathervectors(u);
    s = size(u);
    A = zeros(s);
    for l=1:s(3)
        for p=1:s(4)
            d=u.alpha.*u.F{1,3}(l,:).*u.F{1,4}(p,:);
            A(:,:,l,p) = u.F{1,1}*(spdiags(d(:),0,mm,mm))*u.F{1,2}';
        end
    end
elseif u.dim==5
    u = gathervectors(u);
    s = size(u);
    A = zeros(s);
    for p1=1:s(3)
        for p2=1:s(4)
            for p3=1:s(5)
                A(:,:,p1,p2,p3) = u.F{1,1}*(diag(u.alpha.*u.F{1,3}(p1,:).*u.F{1,4}(p2,:).*u.F{1,5}(p3,:)))*u.F{1,2}';
            end
        end
    end
elseif u.dim==6
    u = gathervectors(u);
    s = size(u);
    A = zeros(s);
    for p1=1:s(3)
        for p2=1:s(4)
            for p3=1:s(5)
                for p4=1:s(6)
                    A(:,:,p1,p2,p3,p4) = u.F{1,1}*(diag(u.alpha.*u.F{1,3}(p1,:).*u.F{1,4}(p2,:).*u.F{1,5}(p3,:).*u.F{1,6}(p4,:)))*u.F{1,2}';
                end
            end
        end
    end
elseif u.dim==7
    u = gathervectors(u);
    s = size(u);
    A = zeros(s);
    for p1=1:s(3)
        for p2=1:s(4)
            for p3=1:s(5)
                for p4=1:s(6)
                    for p5=1:s(7)
                        A(:,:,p1,p2,p3,p4,p5) = u.F{1,1}*(diag(u.alpha.*u.F{1,3}(p1,:).*u.F{1,4}(p2,:).*u.F{1,5}(p3,:).*u.F{1,6}(p4,:).*u.F{1,7}(p5,:)))*u.F{1,2}';
                    end
                end
            end
        end
    end
else
    error('pas programme')
end
