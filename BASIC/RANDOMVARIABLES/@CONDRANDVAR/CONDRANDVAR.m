function rv = CONDRANDVAR(a,varargin)
% function rv = CONDRANDVAR(X,X1,X2,...,Xn,param1,param2)
% X : function_handle pointant sur une RANDVAR
% X1, X2 ..., Xn : RANDVAR auxquelles est conditionnee la variable X
%            les RANDVAR  doivent etre numerotees pour designer une
%            dimension stochastique (utiliser setnumber(Xi,i))
% param1, param2 ... : function_handle parametres de la variable X
% exemple : 
%     CONDRANDVAR(@RVNORMAL,X1,X2,@(X1,X2) X1+X2,@(X1,X2) X1.*X2)
%     variable normale dont la moyenne est X1+X2 et l'ecart-type X1*X2

if ~isa(a,'function_handle')
    error(['first argument must be function_handle'])
end
rv.X = a ;

rv.Y = cell(0,1);
rv.funparam=cell(0,1);
vararg = cell(1,0);
for i=1:length(varargin)    
    if isa(varargin{i},'RANDVARS')
        error('arguments RANDVARS pas acceptes')
    elseif isa(varargin{i},'RANDVAR') || isa(varargin{i},'CONDRANDVAR')
        rv.Y = [rv.Y , varargin(i)];
        vararg = [vararg , {inputname(i+1)}];
    elseif isa(varargin{i},'function_handle') || isa(varargin{i},'inline')
        rv.funparam = [rv.funparam , {fcnchk(varargin{i})}];    
    else
        rv.funparam = [rv.funparam , {(varargin{i})}];        
    end
end



if isa(rv.Y{1},'CONDRANDVAR')
    error('la premiere variable ne doit pas etre conditionnelle')
end

for i=1:length(rv.funparam)
    if isa(rv.funparam{i},'inline') || isa(rv.funparam{i},'function_handle')

        arg = symvar(func2str(rv.funparam{i}));   
        if ~all(ischarin(vararg,arg))
            error(['les fonctions associees aux parametres doivent avoir ' num2str(length(rv.Y)) ' arguments = nombre de variables aleatoires '])    
        end
    end
end

rv.number=[];
rv = class(rv,'CONDRANDVAR');
superiorto('RANDVARS','RANDVAR')

