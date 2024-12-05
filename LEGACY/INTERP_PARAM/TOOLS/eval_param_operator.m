function Ak=eval_param_operator(A,k)
%  function Ak=eval_param_operator(A,k)
%
% A : LRTENSOR with CANONICAL_CORE + TSPACE_OPERATORS
% Need : A.order=2
%        A.space.u{2} : diagonal operator

if length(k)>1
    Ak=cell(length(k),1);
    for i=1:length(k)
        Ak{i}=eval_param_operator(A,k(i));
    end
else
    Ak = A.space.u{1}{1} * A.space.u{2}{1}(k,k);
    for r=1:A.rank
        Ak = Ak + A.space.u{1}{r} * A.space.u{2}{r}(k,k);
    end
end



