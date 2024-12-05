function u = seprand(n,m)
% function u = seprand(n,m)

% global cas

if nargin==1
    m=1;
end
u.dim = length(n);
u.m = m;
u.alpha = ones(1,u.m);
u.F = cell(u.m,u.dim);

for i=1:u.m
    for k=1:u.dim
%u.F{i,k} = (1+rand([n(k),1])).*sign(rand([n(k),1])-1/2);     
% u.F{i,k} = ones(n(k),1);%(1+.1001*rand([n(k),1]));  
%u.F{i,k} = (rand([n(k),1]))
% u.F{i,k} = u.F{i,k}/norm(u.F{i,k});
        u.F{i,k} = randn([n(k),1]);
        u.F{i,k} = u.F{i,k}/norm(u.F{i,k});
    end
end

u = SEPMATRIX(u);

