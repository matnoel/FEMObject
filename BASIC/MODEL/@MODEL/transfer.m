function q2 = transfer(S1,S2,q1,varargin)
% function q2 = transfer(S1,S2,q1)
% transfere la solution q1 definie pour le maillage S1 sur le maillage S2
%
% function q2 = transfer(S1,S2,q1,'nbnodemax',nmax)
% nmax = nombre de noeud maxi traite simultanement (500 par defaut)

if getnbddl(S2)==0
    error('le model n''a pas de ddl -> finaliser le modele')
end

if isa(q1,'PCRADIALMATRIX')
    
    V2 = transfer(S1,S2,double(getV(q1)),varargin{:});
    
    q2 = PCRADIALMATRIX(V2,[S2.nbddl,1],getL(q1));
    return
end

node2 = S2.node;
group2 = getgroupddl(node2);
if isa(q1,'double')
    q2 = zeros(S2.nbddl,size(q1,2));
elseif isa(q1,'PCMATRIX')
    q2 = zeros(S2.nbddl,size(q1,2),getPC(q1));
end

nmax = getcharin('nbnodemax',varargin,500);
for k=1:length(group2)
    subgroup = splitgroup(group2{k},nmax);
    for kk=1:length(subgroup)
        q2temp = eval_sol(S1,q1,POINT(node2(subgroup{kk}.node)),subgroup{kk}.ddlnode);
        numddl = findddl(S2,subgroup{kk}.ddlnode,double(subgroup{kk}.node));
        if isa(q2temp,'double')
            q2temp = double(q2temp);
            q2temp = reshape(q2temp,size(q2,2),length(numddl))';
            q2(numddl,:) = q2temp;
        end
        q2(numddl,:) = q2temp;
    end
end


function G = splitgroup(g,n)
% function G = splitgroup(g,n)

nbnode = size(g.node,1);
ng = floor(nbnode/n);

for k=1:ng
    G{k} = g;
    G{k}.node = g.node((k-1)*n+1:k*n,:);
    G{k}.numddl = g.numddl((k-1)*n+1:k*n,:);
end

if nbnode>n*ng
    k = ng+1;
    G{k} = g;
    G{k}.node = g.node(n*ng+1:end,:);
    G{k}.numddl = g.numddl(n*ng+1:end,:);
end

return
