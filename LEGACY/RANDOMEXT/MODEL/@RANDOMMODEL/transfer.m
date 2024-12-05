function q2 = transfer(S1,S2,q1)
% function q2 = transfer(S1,S2,q1)
% transfere la solution q1 definie pour le maillage S1 sur le maillage S2

if getnbddl(S2)==0
    error('le model n''a pas de ddl -> finaliser le modele')
end

node2 = S2.node;
group2 = getgroupddl(node2);
if isa(q1,'double')
    q2 = zeros(S2.nbddl,size(q1,2));
elseif isa(q1,'PCMATRIX')
    q2 = zeros(S2.nbddl,size(q1,2),getPC(q1));
end

for k=1:length(group2)
    q2temp = eval_sol(S1,q1,POINT(node2(group2{k}.node)),group2{k}.ddlnode);
    numddl = findddl(S2,group2{k}.ddlnode,double(group2{k}.node));
    if isa(q2temp,'double')
        q2temp = double(q2temp);
        q2temp = q2temp(:);
    end
    q2(numddl,:) = q2temp;
end

