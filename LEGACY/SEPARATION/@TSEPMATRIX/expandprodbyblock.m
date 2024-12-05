function ue = expandprodbyblock(A,u,varargin)
% function ue = expandprodbyblock(A,u)
% function ue = expandprodbyblock(A,u,maxblocksize)
% A: TSEPMATRIX
% u: TSEPMATRIX
%
% format of maxblocksize:
%       maxblocksize={10,12,13};
%
%       If maxblocksize={10,0,5}, then maxblocksize{2}=u.m(2);
%
% if nargin==1, maxblocksize={10,10,10}



dim=u.dim;
if nargin==2
    maxblocksize=num2cell(max(u.m,10*ones(1,dim)));
else
    maxblocksize=varargin{1};
end

zerosz=cellfun(@(x) x == 0,maxblocksize);
if any(zerosz)
    maxblocksize{zerosz}=u.m(zerosz);
end

for i=1:dim
    maxblocksize{i}=1:maxblocksize{i}:u.m(i);
    if maxblocksize{i}(end) ~= u.m(i)
        maxblocksize{i} = [maxblocksize{i} u.m(i)];
    end
end
onesblock=cellfun(@(x) numel(x) == 1,maxblocksize);
if any(onesblock)
    maxblocksize(onesblock)={[1 1]};
end

szue=zeros(1,dim);
for d=1:dim
    szue(d)=size(A.F{d}{1},1);
end
ue=zeros(szue);

switch dim
    case 2
        for i=1:numel(maxblocksize{1})-1
            for j=1:numel(maxblocksize{2})-1
                modes1=getmode(maxblocksize{1},i);
                modes2=getmode(maxblocksize{2},j);
                modes={modes1,modes2};
                ue = ue + expand(A*truncate(u,modes));
            end
        end
    case 3
        for i=1:numel(maxblocksize{1})-1
            for j=1:numel(maxblocksize{2})-1
                for k=1:numel(maxblocksize{3})-1
                    modes1=getmode(maxblocksize{1},i);
                    modes2=getmode(maxblocksize{2},j);
                    modes3=getmode(maxblocksize{3},k);
                    modes={modes1,modes2,modes3};
                    ue = ue + expand(A*truncate(u,modes));
                end
            end
        end
    otherwise
        error('not implemented')
end

function mode = getmode(maxblocksz,ind)
    if ind == 1
        mode=maxblocksz(ind):maxblocksz(ind+1);
    else
        mode=maxblocksz(ind)+1:maxblocksz(ind+1);
    end
return
