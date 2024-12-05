function  ae = subsref(a,s)
% function  ae = subsref(a,s)

if length(s)~=2 || ~(isa(s(1).subs{1},'MODEL') || isa(s(1).subs{1},'SEPMODEL'))
    error(['syntaxe : ' inputname(1) '{S,args}(u,...) avec S un (SEP)MODEL et args des arguments supplementaire'])
end

args = s(2).subs;
if length(args)<a.n
    error('rentrer le bon nombre d''arguments')
end

repchar = [];
for j=1:a.n
    if isa(args{j},'char')
        repchar = [repchar,j];
        args{j} = [];
    end
end



if isempty(repchar)
    ae = eval(a,s(1).subs{1},args{:},s(1).subs{2:end});
elseif length(repchar)==1
    ae = calc_vector(a,s(1).subs{1},args{:},s(1).subs{2:end});
elseif length(repchar)==2
    ae = calc_matrix(a,s(1).subs{1},args{:},s(1).subs{2:end});
elseif length(repchar)==3
    ae = calc_trimatrix(a,s(1).subs{1},args{:},s(1).subs{2:end});
else
    error('pas prevu')
end
