function ue = expandbyblock(u,varargin)
% function ue = expandbyblock(u)
% function ue = expandbyblock(u,maxblocksize)
%
% format of maxblocksize:
%       maxblocksize={100,120,130};
%
%       If maxblocksize={10,0,5}, then maxblocksize{2}=u.m(2);
%
% if nargin==1, maxblocksize={100,100,100}

dim=u.dim;

if nargin==1
    maxblocksize=num2cell(100*ones(1,dim));
else
    maxblocksize=varargin{1};
end

zerosz=cellfun(@(x) x == 0,maxblocksize);
if any(zerosz==1)
    maxblocksize{zerosz}=u.m(zerosz);
end

for i=1:dim
    maxblocksize{i}=1:maxblocksize{i}:u.m(i);
    if maxblocksize{i}(end) ~= u.m(i)
        maxblocksize{i} = [maxblocksize{i} u.m(i)];
    end
end

ue=zeros(size(u));
switch dim
    case 2
        for i=1:numel(maxblocksize{1})-1
            for j=1:numel(maxblocksize{2})-1
                modes1=getmode(maxblocksize{1},i);
                modes2=getmode(maxblocksize{2},j);
                modes={modes1,modes2};
                ue = ue + expand(truncate(u,modes));
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
                    ue = ue + expand(truncate(u,modes));
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

