function an = norm(a,varargin)
% function an = norm(a,varargin)

if nargin==1 || (nargin==2 && ~istime(varargin{1})) || (nargin==3 && ~istime(varargin{1}) && ~istime(varargin{2}))
    a = orthogonalizeL(a);
    vn = double(prodscal(a.V,a.V,varargin{:}));
    vn = reshape(vn,a.m,1);
    ln = integratetimes(a.L,a.L);
    an = full(abs(sum(sum(ln.*vn))));
else
    an = full(abs(prodscal(a,a,varargin{:})));
end

if an<100*eps
    an = an;
else
    an = sqrt(an);
end

% if nargin==2 & israndom(varargin{1})
%     an = sqrt(prodscal(a,varargin{1}*a));
% else
%     an = norm(expand(a),varargin{:});
%     try
%         if nargin==2
%             an = full(sqrt(sum(sum(expect(a.*(varargin{1}*a))))));
%         else
%             an = full(sqrt(sum(sum(expect(a.*a)))));
%         end
%         warning('utilisation de la norme PCRADIAL')
%     catch
%         an=norm(expand(a),varargin{:})   ;
%     end
% end
