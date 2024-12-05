function w = mtimes(u,v)

if isa(u,'double')
    if length(u)==1
        w=v;
        w.alpha = v.alpha*u;
    else
        w=v;
        for i=dim
            for j=v.m(i)
                w.F{i}{j}=u*v.F{i}{j};
            end
        end
    end
elseif isa(v,'double')
    w=v*u;
else
    w=u;
    w.m = u.m.*v.m;
    if all(u.m==1)
        w.alpha=double(u.alpha)*v.alpha;
    elseif all(v.m==1)
        w.alpha=double(v.alpha)*u.alpha;
    else
        if w.dim==1
            w.alpha=tensor(kron(double(v.alpha),double(u.alpha)),w.m);
        elseif w.dim==2
            w.alpha=tensor(kron(double(v.alpha),double(u.alpha)),w.m);
        else
            wa=zeros(w.m);
            ua=double(u.alpha);
            va=double(v.alpha);
            vm=v.m;
            um=u.m;
            % matlab can't understand w.alpha([1 2 3])
            % so the algorithm depends on the dimension...
            switch w.dim
                case 3
                    for k=1:vm(3)
                        for j=1:vm(2)
                            for i=1:vm(1)
                                wa(((i-1)*um(1)+1):i*um(1),...
                                    ((j-1)*um(2)+1):j*um(2),...
                                    ((k-1)*um(3)+1):k*um(3))...
                                    =va(i,j,k)*ua;
                            end
                        end
                    end
                    w.alpha=tensor(wa);
                case 4
                    for l=1:v.m(4)
                        for k=1:v.m(3)
                            for j=1:v.m(2)
                                for i=1:v.m(1)
                                    wa(((i-1)*u.m(1)+1):i*u.m(1),...
                                        ((j-1)*u.m(2)+1):j*u.m(2),...
                                        ((k-1)*u.m(3)+1):k*u.m(3),...
                                        ((l-1)*u.m(4)+1):l*u.m(4))...
                                        =v.alpha(i,j,k,l)*u.alpha;
                                end
                            end
                        end
                    end
                    w.alpha=tensor(wa);
                case 5
                    for p=1:v.m(5)
                        for l=1:v.m(4)
                            for k=1:v.m(3)
                                for j=1:v.m(2)
                                    for i=1:v.m(1)
                                        wa(((i-1)*u.m(1)+1):i*u.m(1),...
                                            ((j-1)*u.m(2)+1):j*u.m(2),...
                                            ((k-1)*u.m(3)+1):k*u.m(3),...
                                            ((l-1)*u.m(4)+1):l*u.m(4),...
                                            ((p-1)*u.m(5)+1):p*u.m(5))...
                                            =v.alpha(i,j,k,l,p)*u.alpha;
                                    end
                                end
                            end
                        end
                    end
                    w.alpha=tensor(wa);
                otherwise
                    error('not implemented');
            end
        end
    end

    w.F=cell(w.dim,1);
    for i=1:w.dim
        w.F{i}=cell(w.m(i),1);
    end

    for d=1:w.dim
        c=1;
        for j=1:v.m(d)
            for i=1:u.m(d)
                w.F{d}{c}=u.F{d}{i}*v.F{d}{j};
                c=c+1;
            end
        end
    end
end
