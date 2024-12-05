function [PC,X]=POLYCHAOSSEPDIM(H,p)


if nargin==0

    PC.M = 0;
    PC.p = zeros(0,1);
    PC.n = zeros(0,1);
    PC.PC = cell(1,0);
    h = RANDPOLYS();
    PC = class(PC,'POLYCHAOSSEPDIM',h);
    superiorto('RANDPOLYS')

elseif isa(H,'RANDVARS')
    R = H;
    H = RANDPOLYS(H);
    PC = POLYCHAOSSEPDIM(H,p);
    X = cell(1,getM(R));
    for i=1:getM(R)
        X{i} = project(R{i},PC.PC{i});
    end
    superiorto('RANDPOLYS')
elseif isa(H,'RANDPOLYS')

    PC.M = getM(H);

    if numel(p)==1
        p = repmat(p,1,PC.M);
    end
    for i=1:PC.M
        PC.p(i) = p(i);

        PC.PC{i} = POLYCHAOS(H{i},PC.p(i),'typebase',2);
        PC.n(i) = length(PC.PC{i});
    end

    PC = class(PC,'POLYCHAOSSEPDIM',H);
    superiorto('RANDPOLYS')

    if nargout==2
        R = RANDVARS(H);
        X = cell(1,getM(R));
        for i=1:getM(R)
            X{i} = project(R{i},PC.PC{i});
        end    
    end
end



