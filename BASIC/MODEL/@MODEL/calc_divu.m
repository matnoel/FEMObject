function s = calc_divu(S,q,varargin)
% function s = calc_divu(S,q)
% S : MODEL
% q : solution

if ~isa(S,'MODEL')
    s = calc_divu(q,S,varargin{:});
    return
end

qtemp = q;
if isa(q,'TIMEMATRIX')
    timeq = 1;
    T = gettimemodel(q);
    q = getvalue(q);
elseif isa(q,'TIMERADIALMATRIX')
    timeq = 1;
    T = gettimemodel(q);
    q = double(getV(q));
else
    timeq = 0;
end

if israndom(S) || israndom(q)
    s = calc_divu(S,q,varargin{:});
else
    q = unfreevector(S,q);
    if ischarin('node',varargin)
        s = calc_elemfield(S,@divunode,q,varargin{:});
    else
        if isa(S,'LSMODEL') && getnblevelsets(S)>0
            s = calc_elemfield(S,@divuls,getlevelsets(S),q,varargin{:});
        else
            s = calc_elemfield(S,@divu,q,varargin{:});
        end
        
    end
end

if isa(qtemp,'TIMEMATRIX')
    s = TIMEMATRIX(s,T);
elseif isa(qtemp,'TIMERADIALMATRIX')
    s = TIMERADIALMATRIX(s,T,getL(qtemp));
end
