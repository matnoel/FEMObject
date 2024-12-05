function PC=calc_p_n(PC)
% function PC=calc_p_n(PC)

h = PC.RANDPOLYS;

if isempty(PC.p)
    PC.p=1;
end
if length(PC.p)==1
    PC.p=repmat(PC.p,1,max(1,PC.M));
end

PC.n = PC.p+1;
for k=1:PC.M
    if isa(h{k},'POLYFE')
        PC.n(k) = PC.n(k)*getparam(h{k},'n');
    elseif isa(h{k},'POLYLAGRANGE')
        PC.p(k) = getnbpoints(h{k})-1;
        PC.n(k) = p(k)+1;
    elseif isa(h{k},'POLYFELAGRANGE')
        PC.p(k) = getparam(h{k},'m')-1;
        PC.n(k) = getparam(h{k},'n')*getparam(h{k},'m')-(getparam(h{k},'n')-1);
    elseif isa(h{k},'POLYWAVELETS')
        if getparam(h{k},'p')~=PC.p(k)
            h{k} = setnumber(POLYWAVELETS(getparam(h{k},'n'),PC.p(k)),getnumber(h{k}));
        end
        PC.n(k) = 2^getparam(h{k},'n')*PC.n(k);
    end
end

% for k=1:PC.M
%     PC.n(k) = length(unique(PC.indices(:,k)));
% end
