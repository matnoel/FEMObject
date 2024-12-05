function [u,flag,result,bini] = dsolve_reuse(GSD,T,A,B,b,ureuse,varargin)
% function [u,flagreuse,result,bini] = dsolve_reuse(GSD,T,A,B,b,ureuse,varargin)
% resolution de Adot(u)+Bu=b sur la base reduite issue de ureuse
% attention, la prise en compte de la condition ini doit etre dans b
% ureuse : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)
% u : solution : PCRADIALMATRIX 
% bini = b-A*dot(u)-Bu (residu)


paramradial.reuse = getparam(GSD,'reuse');
if paramradial.reuse && ~isempty(ureuse)
   if isa(ureuse,'PCTIMEMATRIX') && ispcradial(ureuse)
        if length(ureuse)==0
            paramradial.reuse=false;
        end
   else
       fprintf('  utiliser une PCTIMEMATRIX sous forme radial ou une MULTIMATRIX\n')
            paramradial.reuse=false;
   end
else
    paramradial.reuse = false;
end

n=size(A,1);
nt = length(T);

if paramradial.reuse

   if isa(ureuse,'PCTIMEMATRIX') && ispcradial(ureuse)
            ureuse = getV(ureuse);
   end
   
toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');

if directsolve
localstosolver = @(aU,fU) solve(expand(aU),fU);
else
localstosolver = @(aU,fU) cgs(aU,fU,toliter,[],'noupdate');
end


ureuse = getvalue(ureuse);
if isa(ureuse,'MULTIMATRIX')
    V = cell(1,length(ureuse));
    DV = cell(1,length(ureuse));
   for i=1:length(ureuse)
    V{i} = TIMEMATRIX(ureuse{i},T,[n,1]);   
    DV{i} = diff(V{i});
   end
elseif isa(ureuse,'double')
    V ={ureuse};
end
       
       
m = length(V);
fU = cell(1,m);
aU = cell(m,m);
bU = cell(m,m);

for i=1:length(V)
fU{i} = integratemtimes(V{i}',b);
for j=1:length(V)
aU{i}{j} = integratemtimes(V{i}',A,DV{j});
bU{i}{j} = integratemtimes(V{i}',B,V{j});
aU{i}{j} = aU{i}{j}+bU{i}{j};
end
aU{i} = horzcat(aU{i}{:});
end
fU = vertcat(fU{:});
aU = vertcat(aU{:});

l = localstosolver(aU,fU);

result.Rayg = full(expect(fU,l));
if isa(result.Rayg,'MULTIMATRIX')
result.Rayg = double(cell2mat(result.Rayg));
end

result.rayg = trace(result.Rayg);

bini=b;

for i=1:length(V)
Vi = PCRADIALMATRIX(getvalue(V{i}),[n,nt],l(i));
if i==1
u = PCTIMEMATRIX(Vi,T,[n,1]);    
else
u= u + PCTIMEMATRIX(Vi,T,[n,1]);
end

if nargout>=4 
bini = bini - (A*DV{i}+B*V{i})*l(i);
end

end

else
u = PCTIMEMATRIX(PCRADIALMATRIX([n,nt],getPC(b)),T,[n,1]);
result.Rayg = 0;
result.rayg = 0;
bini = b;
end

flag= ~paramradial.reuse;

