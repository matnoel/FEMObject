function lsx=eval(ls,x,varargin)

for i=1:length(varargin)
    varargin{i}=MYDOUBLEND(reshape(full(varargin{i}),[1,1,1,numel(varargin{i})]));
end

P1 = POINT([varargin{1:2}]);
P2 = POINT([varargin{3:4}]);
c = (P1+P2)/2;
r = (P2-P1)/2;
lsx = abs(x - c) - r;

lsx = max(lsx',[],2);

