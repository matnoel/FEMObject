function varargout = dilateaxis(fact,varargin)
if nargin==1
    ax=axis();
else
    ax = varargin{1};   
end
n = length(ax)/2;
if numel(fact)==1
    fact = repmat(fact,1,n);
end
for i=1:n
    am = (ax(2*i-1)+ax(2*i))/2;
    da = (ax(2*i)-ax(2*i-1));
    ax(2*i-1:2*i) = [am-da/2*fact(i),am+da/2*fact(i)];
end

if nargout==1
    varargout{1}=ax;
else
    axis(ax);
end
